module ParamControl

STATUS_ACCEPT   = 0
STATUS_CONTINUE = 1
STATUS_FAIL     = 2


mutable struct param_controller{T <: AbstractFloat}
    beg_i           :: Integer
    end_i           :: Integer
    params          :: Array{T}
    test_param      :: Array{T}
    saved_param     :: Array{T}
    fail_count      :: Integer
    fail_count_max  :: Integer
    N               :: Integer

    function param_controller(
        params::Array{T},
        fail_count_max::Integer
    ) where T <: AbstractFloat

        return new{T}(
            1,
            1,
            params,
            copy(params[1,:]),
            copy(params[1,:]),
            0,
            fail_count_max,
            size(params)[1],
        )
    end
end

function iterate_and_adjust!(
    o::param_controller,
    converge_flag::Bool
)

    if converge_flag == true
        o.fail_count = 0
        if o.test_param == o.params[o.end_i, :]
            #println("Touch endpoint")
            # If now touches the end point
            
            if o.end_i == o.N
                # If it is entirely complete
                return STATUS_ACCEPT
            else
                # If not then move on to next interval
                o.beg_i, o.end_i = o.end_i, o.end_i + 1
                o.test_param[:] = o.params[o.end_i, :]
                return STATUS_CONTINUE
            end
        else
            #println("Does not touch endpoint")
            #println(o.ab_pair)
            # If not touching end point
            o.saved_param[:] = o.test_param
            o.test_param[:]  = o.params[o.end_i, :]

            return STATUS_CONTINUE
        end
    else
        # If not converge
        o.fail_count += 1
        if o.fail_count >= o.fail_count_max
            return STATUS_FAIL
        else
            # Backup for a little
            o.test_param[:] = (o.saved_param[:] + o.test_param) / 2.0
            #println("Fail count: ", fail_count, "; ab_pair change to :", ab_pair)
            return STATUS_CONTINUE
        end
    end
end


end 
