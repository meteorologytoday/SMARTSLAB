function OC_doConvectiveAdjustment!(oc::OceanColumn)
    if_adjust, oc.b_ML, oc.h_ML, oc.FLDO = doConvectiveAdjustment!(
        zs   = oc.zs,
        bs   = oc.bs,
        h_ML = oc.h_ML,
        b_ML = oc.b_ML,
        FLDO = oc.FLDO,
    )

    return if_adjust
end


"""

This function only do convective adjustment for the upper most mixed layer.
It searches for the lowest layer that has larger buoyancy than mixed-layer then mixed all layers above it.

"""
function doConvectiveAdjustment!(;
    zs   :: Array{Float64, 1},
    bs   :: Array{Float64, 1},
    h_ML :: Float64,
    b_ML :: Float64,
    FLDO :: Integer,
)
    
    if_adjust = false

    if FLDO == -1
        return if_adjust, b_ML, h_ML, FLDO 
    end

    # 1. Search from bottom to see if buoyancy is monotically increasing
    # 2. If not, record the peak, then keep detecting until hitting the 
    #    layer having equal or greater buoyancy. Record this interval.
    # 3. Find the minimum value in this interval b_min.
    # 4. Use b_min to decide the bottom layer going to be mixed.
    # 5. Mix this interval.

    new_b_ML = b_ML
    new_h_ML = h_ML
    new_FLDO = FLDO

    stage = :reset
    peak_layer = 0
    top_layer = 0
    bot_layer = 0
    b_peak = 0.0



   for i = length(bs):-1:FLDO


        if stage == :reset
            peak_layer = 0
            top_layer = 0
            bot_layer = 0
            b_peak = 0.0
            stage = :search_peak_layer
        end

        if stage == :search_peak_layer

            Δb = bs[i] - ((i==FLDO) ? b_ML : bs[i-1])
            #println("FLDO:", FLDO, "; i:", i, "; Δb:", Δb)
            if Δb > 0.0
                if_adjust = true
                stage = :search_top_layer
                peak_layer = i
                b_peak = bs[peak_layer]
            else
                continue
            end
        end

        if stage == :search_top_layer

            #println(":search_top_layer")
            if i == FLDO
                top_layer = (b_ML > b_peak) ? FLDO : -1
                stage = :search_bot_layer
            elseif bs[i-1] > b_peak
                top_layer = i
                stage = :search_bot_layer
            else
                continue
            end
        end

        if stage == :search_bot_layer

            #println(":search_bot_layer")

            if peak_layer == length(bs)

                bot_layer = peak_layer
                stage = :start_adjustment

            else
                b_min = 0.0
                if top_layer == -1
                    b_min = min(b_ML, minimum(bs[FLDO:peak_layer]))
                else
                    b_min = minimum(bs[top_layer:peak_layer])
                end 

                bot_layer = peak_layer + 1
                while true
                    if bs[bot_layer] >= b_min
                        if bot_layer == length(bs)
                            stage = :start_adjustment
                            break
                        else
                            bot_layer += 1
                        end
                    else
                        stage = :start_adjustment
                        break
                    end
                end
            end 
        end


        if stage == :start_adjustment
            #println(":start_adjustment")
            bot_z = zs[bot_layer+1]
            top_z = (top_layer == -1) ? 0.0 : (
                 (top_layer == FLDO) ? -h_ML : zs[top_layer]
            )

            new_b = (getIntegratedBuoyancy(
                zs       =  zs,
                bs       =  bs,
                b_ML     =  b_ML,
                h_ML     =  h_ML,
                target_z =  bot_z
            ) - getIntegratedBuoyancy(
                zs       =  zs,
                bs       =  bs,
                b_ML     =  b_ML,
                h_ML     =  h_ML,
                target_z =  top_z
            ))  / (top_z - bot_z)
           
        
            if top_layer == -1
                new_b_ML = new_b
                new_h_ML = (bot_layer == length(bs)) ? zs[1] - zs[end] : zs[1] - zs[bot_layer + 1]
                new_FLDO = setMixedLayer!(bs=bs, zs=zs, b_ML=new_b_ML,h_ML=new_h_ML)
            else
                bs[top_layer:bot_layer] .= new_b
            end 

            stage = :reset 
        end

    end

    return if_adjust, new_b_ML, new_h_ML, new_FLDO
end

           

