import jmespath.parser
import jmespath.visitor
import jmespath.functions
import itertools
import types

# Register generator type in jmespath
jmespath.functions.TYPES_MAP['generator'] = 'array'
jmespath.functions.REVERSE_TYPES_MAP['array'] += ('generator',)


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


class TreeInterpreterGenerator(jmespath.visitor.TreeInterpreter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._generators = {}

    def _gen_to_list(self, gen):
        """
        Called when a jmespath operation requires random access to generator.
        Memoises already converted generators and returns the already generated list.
        :param gen: generator object
        :return: list of values returned by generator
        """
        gen = self._generators.get(id(gen), gen)
        if isinstance(gen, types.GeneratorType):
            l = list(gen)
            self._generators[id(gen)] = l
            return l
        return gen

    def visit(self, node, *args, **kwargs):
        # if a visit caused list conversion, get list. Assume that 'value' is args[0].
        if len(args) and isinstance(args[0], types.GeneratorType):
            args = list(args) #convert from tuple
            args[0] = self._generators.get(id(args[0]), args[0])
        return super().visit(node, *args, **kwargs)

    def visit_field(self, node, value):
        try:
            return value.get(node['value'])
        except AttributeError:
            # Allow accessing objects fields TODO push this change upstream, possibly with config flag
            try:
                return getattr(value, node['value'])
            except AttributeError:
                return None

    def visit_function_expression(self, node, value):
        resolved_args = []
        for child in node['children']:
            current = self._gen_to_list(self.visit(child, value))
            resolved_args.append(current)
        return self._functions.call_function(node['value'], resolved_args)

    def visit_not_expression(self, node, value):
        original_result = self.visit(node['children'][0], value)
        if original_result == 0:
            # Special case for 0, !0 should be false, not true.
            # 0 is not a special cased integer in jmespath.
            return False
        return self._is_false(original_result) # TODO bugfix, push this change upstream

    def visit_filter_projection(self, node, value):
        base = self.visit(node['children'][0], value)
        if not (isinstance(base, list) or isinstance(base, types.GeneratorType)):
            return None
        comparator_node = node['children'][2]
        for element in base:
            comparison = self.visit(comparator_node, element)
            if self._is_true(comparison):
                current = self.visit(node['children'][1], element)
                if current is not None:
                    yield current

    def visit_flatten(self, node, value):
        base = self.visit(node['children'][0], value)
        if not isinstance(base, list):
            # Can't flatten the object if it's not a list.
            return None
        for element in base:
            if isinstance(element, list):
                for subelement in element:
                    yield subelement
            else:
                yield element

    def visit_index(self, node, value):
        value = self._gen_to_list(value)
        return super().visit_index(node, value)

    def visit_slice(self, node, value):
        return itertools.islice(value, *node['children'])

    def visit_multi_select_list(self, node, value):
        if value is None:
            return None
        for child in node['children']:
            yield self.visit(child, value)

    def visit_projection(self, node, value):
        base = self.visit(node['children'][0], value)
        if not isinstance(base, list):
            return None
        for element in base:
            current = self.visit(node['children'][1], element)
            if current is not None:
                yield current

    def visit_value_projection(self, node, value):
        base = self.visit(node['children'][0], value)
        try:
            base = base.values()
        except AttributeError:
            return None
        for element in base:
            current = self.visit(node['children'][1], element)
            if current is not None:
                yield current

    def _is_false(self, value):
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
