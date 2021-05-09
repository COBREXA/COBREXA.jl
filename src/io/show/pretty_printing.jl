"""
    _print_with_colors(io, def::String, payload; kwargs...)

Prints nicely colorized keys and values.
"""
_print_with_colors(io, def::String, payload; kwargs...) =
    _print_with_colors(io, def, isnothing(payload) ? "---" : string(payload); kwargs...)

"""
    _print_with_colors(
        io,
        def::String,
        payload::String;
        def_color = _constants.colors.key,
        empty_color = _constants.colors.empty,
        payload_color = _constants.colors.payload,
    )

Specialization of `_print_with_colors` for plain strings.
"""
function _print_with_colors(
    io,
    def::String,
    payload::String;
    def_color = _constants.colors.key,
    empty_color = _constants.colors.empty,
    payload_color = _constants.colors.payload,
)
    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = empty_color), "---")
    else
        println(io, Crayon(foreground = payload_color), payload)
    end
end

"""
    _print_with_colors(
        io,
        def::String,
        payload::Dict;
        def_color = _constants.colors.key,
        empty_color = _constants.colors.empty,
        payload_color = _constants.colors.payload,
    )

Specialization of `_print_with_colors` for dictionaries.
"""
function _print_with_colors(
    io,
    def::String,
    payload::Dict;
    def_color = _constants.colors.key,
    empty_color = _constants.colors.empty,
    payload_color = _constants.colors.payload,
)

    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = empty_color), "---")
    else
        println(io, "")
        for (k, v) in payload
            if length(v) > 2 && length(v[1]) < 20
                println(
                    io,
                    Crayon(foreground = payload_color),
                    "\t",
                    k,
                    ": ",
                    v[1],
                    ", ..., ",
                    v[end],
                )
            elseif length(v[1]) > 20 # basically for envipath annotations... or long notes
                println(
                    io,
                    Crayon(foreground = payload_color),
                    "\t",
                    k,
                    ": ",
                    v[1][1:20],
                    "...",
                )
            else
                println(io, Crayon(foreground = payload_color), "\t", k, ": ", v)
            end
        end
    end
end