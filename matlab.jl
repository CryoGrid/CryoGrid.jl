module matlab
    using Interpolations

    function interp1(xpt, ypt, x, method)
        extrapvalue = nothing
        if extrapvalue == nothing
            y = zeros(x)
            idx = trues(x)
        else
            y = extrapvalue*ones(x)
            idx = (x .>= xpt[1]) .& (x .<= xpt[end])
        end

        if method == "linear"
            intf = interpolate((xpt,), ypt, Gridded(Linear()))
            y[idx] = intf[x[idx]]

        elseif method == "cubic"
            itp = interpolate(ypt, BSpline(Cubic(Natural())), OnGrid())
            intf = scale(itp, xpt)
            y[idx] = [intf[xi] for xi in x[idx]]
        elseif method == "nearest"
            intf = interpolate((xpt,), ypt, Gridded(Constant()))
            y[idx] = intf[x[idx]]
        end
        return y
    end

    function nansum(a::Array{Float64,1})
        out=0.0;
        @inbounds @fastmath for i=1:length(a)
            if !isnan(a[i])
                out=out+a[i]
            end
        end
        return out
    end

    function datenum(d)
        #=
             datenum(d::Dates.DateTime)
        Converts a Julia DateTime to a MATLAB style DateNumber.
        MATLAB represents time as DateNumber, a double precision floating
        point number being the the number of days since January 0, 0000
        Example
            datenum(now())
        =#
        out = Array{Float64,1}(length(d))
        MATLAB_EPOCH = Dates.DateTime(-0001,12,31)
        fac = 1000. * 60. * 60. * 24.
        for i=1:length(d)
            out[i] = Dates.value(d[i] - MATLAB_EPOCH) / fac
        end
        return out
    end

    function datestr(d::Array{Float64,1})
        datestr=Array{DateTime,1}(length(d))
        MATLAB_EPOCH = Dates.DateTime(-0001,12,31)
        fac = 1000. * 60. * 60. * 24.
        for i=1:length(d)
            datestr[i] = Dates.Millisecond(d[i]*fac) + MATLAB_EPOCH
        end
        return datestr
    end

    function fastminmax(a_min,a_max,b)
        #min_out=Array{Float64,1}(length(b))
        #max_out=Array{Float64,1}(length(b))
        @inbounds @fastmath for i=1:length(b)
            if b[i]<a_min[i]
                a_min[i] = b[i]
            else
                a_min[i] = a_min[i]
            end
            if b[i]>a_max[i]
                a_max[i] = b[i]
            else
                a_max[i] = a_max[i]
            end
        end
        return a_min, a_max
    end

end
