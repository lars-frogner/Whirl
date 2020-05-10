const DEBUG = true

"""
    @dbgassert(expr, [message]) -> true

If `DEBUG = true`, throws an error with an optionally specified `message` if the
given code expression evaluates to `false`, otherwise returns `true`.
"""
macro dbgassert(expr, msgs...)
    if DEBUG
        msg_str = isempty(msgs) ? string(expr) : msgs[1]
        esc(:($expr || throw(AssertionError($msg_str))))
    else
        true
    end
end

"""
    @dbgasserteq(left, right, [message]) -> true

If `DEBUG = true`, throws an error with an optionally specified `message` if the
given `left` code expression does not evaluates to the same value as the given
`right` expression, otherwise returns `true`.
"""
macro dbgasserteq(left, right, msgs...)
    if DEBUG
        msg_str = isempty(msgs) ? "$(left) == $(right)" : msgs[1]
        esc(:($left == $right || throw(AssertionError($msg_str))))
    else
        true
    end
end
