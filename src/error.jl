"""
    @check(expr, [message]) -> true

Throws an error with an optionally specified `message` if the given code
expression evaluates to `false`, otherwise returns `true`.
"""
macro check(expr, msgs...)
    msg_str = isempty(msgs) ? string(expr) : msgs[1]
    esc(:($expr || throw(AssertionError($msg_str))))
end

"""
    @checkeq(left, right, [message]) -> true

Throws an error with an optionally specified `message` if the given `left` code
expression does not evaluates to the same value as the given `right`
expression, otherwise returns `true`.
"""
macro checkeq(left, right, msgs...)
    msg_str = isempty(msgs) ? "$(left) == $(right)" : msgs[1]
    esc(:($left == $right || throw(AssertionError($msg_str))))
end
