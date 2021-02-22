# workaround for @test not allowing error messages
printerror(str) = begin (@error str); false; end

function allfinite(x::AbstractArray)
	chk = isfinite.(x)
	all(chk) || printerror("Infinite values at $(findall(.!chk)) in $x")
end

function allequal(x::AbstractArray, y::AbstractArray; atol=1.0e-3)
	chk = .≈(x,y,atol=atol)
	all(chk) || printerror("Values do not match at $(findall(.!chk)) in $x and $y")
end
