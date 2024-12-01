function create_susceptibility_by_year(year::Int, value1::Float64, value2::Float64, value3::Float64)
    # fixed the first four values
    base_susceptibility = [
        fill(value1), #0-4
        fill(1), #5-9
        fill(1), #10-19
        fill(1) #20-29
    ]
    # susceptibility values change depending on the year
    if year == 2010
        append!(base_susceptibility, [fill(value3), fill(value3), fill(value3)])
    elseif year == 2015
        append!(base_susceptibility, [fill(value2), fill(value3), fill(value3)])
    elseif year == 2020
        append!(base_susceptibility, [fill(1), fill(value3), fill(value3)])
    elseif year == 2024
        append!(base_susceptibility, [fill(1), fill(value2), fill(value3)])
    elseif year == 2030
        append!(base_susceptibility, [fill(1), fill(1), fill(value3)])
    else
        error("Invalid year. Provide a year between 2010 and 2030.")
    end
    return base_susceptibility
end

function create_susceptibilities_for_years(years::Vector{Int}, value1::Float64, value2::Float64, value3::Float64)
    susceptibilities = Dict{Int, Vector}()
    for year in years
        susceptibilities[year] = create_susceptibility_by_year(year, value1, value2, value3)
    end
    return susceptibilities
end

# %%
function model_averaging_sexual(s1, s2, s3, s4, dominant_eigvals_all, dominant_eigvals_phys, dominant_eigvals_home, dominant_eigvals_physhome)
    weighted_avg_sexual = []
    w1 = 0.22 #all
    w2 = 0.38 #phys
    w3 = 0.21 #home
    w4 = 0.17 #physhome
    for i in 1:length(dominant_eigvals_all)
        if i > length(dominant_eigvals_all) - 2  
            weighted_value = (1/s1) * w1 * dominant_eigvals_all[i] + (1/s2) * w2 * dominant_eigvals_phys[i] + (1/s3) * w3 * dominant_eigvals_home[i] + (1/s4) * w4 * dominant_eigvals_physhome[i]
        else
            weighted_value = w1 * dominant_eigvals_all[i] + w2 * dominant_eigvals_phys[i] + w3 * dominant_eigvals_home[i] + w4 * dominant_eigvals_physhome[i]
        end
        push!(weighted_avg_sexual, weighted_value)
    end
    return weighted_avg_sexual
end
