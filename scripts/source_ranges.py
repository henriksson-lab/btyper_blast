"""Locate selected unchanged C/C++ definitions; not a dependency-graph parser."""
import re


def definition(text, pattern, start=0, stop=None, semicolon=False):
    match = re.compile(pattern, re.MULTILINE).search(text, start, stop or len(text))
    if not match:
        raise ValueError(f"Missing reference definition: {pattern}")
    opening = text.index("{", match.end())
    depth = 0
    tokens = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}]', re.S)
    for token in tokens.finditer(text, opening):
        if token[0] == "{":
            depth += 1
        elif token[0] == "}":
            depth -= 1
            if depth == 0:
                end = token.end()
                if semicolon:
                    assert text[end] == ";"
                    end += 1
                return match.start(), end
    raise ValueError(f"Unterminated reference definition: {pattern}")
