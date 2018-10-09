struct ηCtl
    N               :: Integer
    η_growth_factor :: AbstractFloat
    η_reset_value   :: AbstractFloat
    η_now           :: AbstractFloat
    ϵ               :: Array{AbstractFloat,1}
    p               :: Integer                 # pointer to the latest ϵ

    ηCtl(η, N) = new(η, N, Inf, 1)
end


function judge(
    ctl :: ηCtl,
    ϵ   :: AbstractFloat
)


end
