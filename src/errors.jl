function notimplemented(::T)
    error("$(StackTraces.stacktrace()[2].func) not implemented for type $(T)")
end
