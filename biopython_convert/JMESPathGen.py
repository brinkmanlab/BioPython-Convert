import jmespath.parser
import jmespath.visitor
import jmespath.functions
import jmespath.exceptions
import itertools
import types

# Register generator type in jmespath
jmespath.functions.TYPES_MAP['generator'] = 'array'
jmespath.functions.REVERSE_TYPES_MAP['array'] += ('generator',)

# Register biopython types in jmespath
jmespath.functions.TYPES_MAP['Seq'] = 'string'
jmespath.functions.REVERSE_TYPES_MAP['string'] += ('Seq',)
jmespath.functions.REVERSE_TYPES_MAP['Seq'] = ('Seq',)
jmespath.functions.REVERSE_TYPES_MAP['SeqFeature'] = ('SeqFeature',)
jmespath.functions.TYPES_MAP['ExactPosition'] = 'number'
jmespath.functions.TYPES_MAP['BeforePosition'] = 'number'
jmespath.functions.TYPES_MAP['BetweenPosition'] = 'number'
jmespath.functions.TYPES_MAP['AfterPosition'] = 'number'
jmespath.functions.TYPES_MAP['OneOfPosition'] = 'number'
jmespath.functions.REVERSE_TYPES_MAP['number'] += ('ExactPosition',)
jmespath.functions.REVERSE_TYPES_MAP['number'] += ('BeforePosition',)
jmespath.functions.REVERSE_TYPES_MAP['number'] += ('BetweenPosition',)
jmespath.functions.REVERSE_TYPES_MAP['number'] += ('AfterPosition',)
jmespath.functions.REVERSE_TYPES_MAP['number'] += ('OneOfPosition',)

# this implementation includes https://github.com/jmespath/jmespath.site/pull/6
# and https://github.com/jmespath/jmespath.py/issues/159


def compile(expression):
    return Parser().parse(expression)


def search(expression, data, options=None):
    return Parser().parse(expression).search(data, options=options)


class Parser(jmespath.parser.Parser):
    def _parse(self, expression):
        result = super()._parse(expression)
        return ParsedResult(result.expression, result.parsed)


class ParsedResult(jmespath.parser.ParsedResult):
    def search(self, value, options=None):
        interpreter = TreeInterpreterGenerator(options)
        result = interpreter.visit(self.parsed, value)
        return result


class ExtendedFunctions(jmespath.functions.Functions):
    def call_function(self, function_name, resolved_args, **kwargs):
        try:
            spec = self.FUNCTION_TABLE[function_name]
        except KeyError:
            raise jmespath.exceptions.UnknownFunctionError(
                "Unknown function: %s()" % function_name)
        function = spec['function']
        signature = spec['signature']
        self._validate_arguments(resolved_args, signature, function_name)
        return function(self, *resolved_args, **kwargs)

    @jmespath.functions.signature({'types': ['object']}, {'types': ['expref']})
    def _func_let(self, lexical_scope, expref, **kwargs):
        if 'scope' in kwargs:
            scope = dict(kwargs['scope'])
            scope.update(lexical_scope)
        else:
            scope = dict(lexical_scope)
        kwargs['scope'] = scope
        return expref.visit(expref.expression, expref.context, **kwargs)

    @jmespath.functions.signature({'types': ['string']}, {'types': ['string']})
    def _func_split(self, on, val):
        return val.split(on)

    @jmespath.functions.signature({'types': ['Seq']}, {'types': ['SeqFeature']})
    def _func_extract(self, seq, feature):
        return feature.extract(seq)


class _Expression(jmespath.visitor._Expression):
    def __init__(self, expression, interpreter, context):
        super().__init__(expression, interpreter)
        self.context = context


class TreeInterpreterGenerator(jmespath.visitor.TreeInterpreter):
    def __init__(self, options=None, *args, **kwargs):
        options = options or jmespath.visitor.Options(custom_functions=ExtendedFunctions())
        super().__init__(*args, options=options, **kwargs)
        self._generators = {}

    def _gen_to_list(self, gen, recurse=False):
        """
        Called when a jmespath operation requires random access to generator.
        Memoises already converted generators and returns the already generated list.
        :param gen: generator object
        :return: list of values returned by generator
        """
        if getattr(gen, '__hash__', None) is not None:
            gen = self._generators.get(gen, gen)
        if isinstance(gen, types.GeneratorType):
            l = list(gen)
            self._generators[gen] = l
            if recurse:
                for i in range(len(l)):
                    if isinstance(l[i], types.GeneratorType):
                        l[i] = self._gen_to_list(l[i])
            return l
        return gen

    def visit(self, node, *args, **kwargs):
        # if a visit caused list conversion, get list. Assume that 'value' is args[0].
        if len(args) and isinstance(args[0], types.GeneratorType):
            args = list(args)  # convert from tuple
            args[0] = self._generators.get(args[0], args[0])
        return super().visit(node, *args, **kwargs)

    def visit_field(self, node, value, **kwargs):
        try:
            return value.get(node['value'])
        except AttributeError:
            # Allow accessing objects fields TODO push this change upstream, possibly with config flag
            try:
                return getattr(value, node['value'])
            except AttributeError:
                # If the field is not defined in the current object, then fall back
                # to checking in the scope chain, if there's any that has been
                # created.
                if 'scope' in kwargs:
                    return kwargs['scope'].get(node['value'], None)
                return None

    def visit_function_expression(self, node, value, **kwargs):
        resolved_args = []
        for child in node['children']:
            current = self._gen_to_list(self.visit(child, value, **kwargs), True)
            resolved_args.append(current)
        return self._functions.call_function(node['value'], resolved_args)

    def visit_not_expression(self, node, value, **kwargs):
        original_result = self.visit(node['children'][0], value, **kwargs)
        if original_result == 0:
            # Special case for 0, !0 should be false, not true.
            # 0 is not a special cased integer in jmespath.
            return False
        return self._is_false(original_result) # TODO bugfix, push this change upstream

    def visit_filter_projection(self, node, value, **kwargs):
        base = self.visit(node['children'][0], value, **kwargs)
        if not isinstance(base, (list, types.GeneratorType, map, filter)):
            return None
        comparator_node = node['children'][2]
        for element in base:
            comparison = self.visit(comparator_node, element, **kwargs)
            if self._is_true(comparison):
                current = self.visit(node['children'][1], element, **kwargs)
                if current is not None:
                    yield current

    def visit_flatten(self, node, value, **kwargs):
        base = self.visit(node['children'][0], value, **kwargs)
        if not isinstance(base, (list, types.GeneratorType, map, filter)):
            # Can't flatten the object if it's not a list.
            return None
        for element in base:
            if isinstance(element, (list, types.GeneratorType, map, filter)):
                for subelement in element:
                    yield subelement
            else:
                yield element

    def visit_index(self, node, value, **kwargs):
        value = self._gen_to_list(value)
        return super().visit_index(node, value)

    def visit_slice(self, node, value, **kwargs):
        return itertools.islice(value, *node['children'])

    def visit_multi_select_list(self, node, value, **kwargs):
        if value is None:
            return None
        for child in node['children']:
            yield self.visit(child, value, **kwargs)

    def visit_projection(self, node, value, **kwargs):
        base = self.visit(node['children'][0], value, **kwargs)
        if not isinstance(base, (list, types.GeneratorType, map, filter)):
            return None
        for element in base:
            current = self.visit(node['children'][1], element, **kwargs)
            if current is not None:
                yield current

    def visit_value_projection(self, node, value, **kwargs):
        base = self.visit(node['children'][0], value, **kwargs)
        try:
            base = base.values()
        except AttributeError:
            return None
        for element in base:
            current = self.visit(node['children'][1], element, **kwargs)
            if current is not None:
                yield current

    def _is_false(self, value, **kwargs):
        if isinstance(value, types.GeneratorType):
            # peek generator instead of _gen_to_list()
            try:
                peek = next(value)
                old_value = value
                value = itertools.chain((peek,), value)
                self._generators[id(old_value)] = value
                return False
            except StopIteration:
                return True
        return super()._is_false(value)

    def visit_expref(self, node, value, **kwargs):
        return _Expression(node['children'][0], self, value)

    def visit_subexpression(self, node, value, **kwargs):
        result = value
        for node in node['children']:
            result = self.visit(node, result, **kwargs)
        return result

    def visit_comparator(self, node, value, **kwargs):
        # Common case: comparator is == or !=
        comparator_func = self.COMPARATOR_FUNC[node['value']]
        if node['value'] in self._EQUALITY_OPS:
            return comparator_func(
                self.visit(node['children'][0], value, **kwargs),
                self.visit(node['children'][1], value, **kwargs)
            )
        else:
            # Ordering operators are only valid for numbers.
            # Evaluating any other type with a comparison operator
            # will yield a None value.
            left = self.visit(node['children'][0], value, **kwargs)
            right = self.visit(node['children'][1], value, **kwargs)
            num_types = (int, float)
            if not (jmespath._is_comparable(left) and
                    jmespath._is_comparable(right)):
                return None
            return comparator_func(left, right)

    def visit_current(self, node, value, **kwargs):
        return super().visit_current(node, value)

    def visit_identity(self, node, value, **kwargs):
        return super().visit_identity(node, value)

    def visit_index_expression(self, node, value, **kwargs):
        result = value
        for node in node['children']:
            result = self.visit(node, result, **kwargs)
        return result

    def visit_key_val_pair(self, node, value, **kwargs):
        return self.visit(node['children'][0], value, **kwargs)

    def visit_literal(self, node, value, **kwargs):
        return super().visit_literal(node, value)

    def visit_multi_select_dict(self, node, value, **kwargs):
        if value is None:
            return None
        collected = self._dict_cls()
        for child in node['children']:
            collected[child['value']] = self.visit(child, value, **kwargs)
        return collected

    def visit_or_expression(self, node, value, **kwargs):
        matched = self.visit(node['children'][0], value, **kwargs)
        if self._is_false(matched):
            matched = self.visit(node['children'][1], value, **kwargs)
        return matched

    def visit_and_expression(self, node, value, **kwargs):
        matched = self.visit(node['children'][0], value, **kwargs)
        if self._is_false(matched):
            return matched
        return self.visit(node['children'][1], value, **kwargs)

    def visit_pipe(self, node, value, **kwargs):
        result = value
        for node in node['children']:
            result = self.visit(node, result, **kwargs)
        return result

    def _is_true(self, value, **kwargs):
        return super()._is_true(value)

