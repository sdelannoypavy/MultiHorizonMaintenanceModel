using NCDatasets
using Dates
using DataFrames
using Statistics
using DelimitedFiles

ds = Dataset("data_merged.nc", "r")

# shww = significant height of wind waves 
shww = ds["shww"][:]

# data as vector
shww_data = [shww[i, 1, 1] for i in 1:size(shww, 1)]

start_date = DateTime("1940-01-01T00:00:00")
n_hours = length(shww_data)
date_range = start_date:Hour(1):(start_date + Hour(n_hours - 1))

# Create dataframe
df = DataFrame(timestamp = date_range, values = shww_data)


# Hourly data -> compute daily mean
df.day = Date.(df.timestamp)
daily_means = combine(groupby(df, :day), :values => mean => :daily_mean)

# create scenarios - extract 2 month of data for each. 
# We create 10 scenarios using 10 years (no data analysis)
scenarios_shww = fill(-1.0, 10, 6, 62)

for S in 1:10
    for T in 1:6
        month_list = [(T-1)*2 + 1, (T-1)*2 + 2]   
        y = 1939 + S 

        # extract data for the corresponding year and strategic period
        subset = filter(row -> year(row.day) == y && month(row.day) in month_list, daily_means)
        shww_by_day_subset = subset.daily_mean

        l = length(shww_by_day_subset)
        scenarios_shww[S,T,1:l] = shww_by_day_subset
    end
end

function Production(shww)

    if shww == -1
        return -1
    else
        wind_speed = -1 + sqrt(1+(0.016/0.406^2)*4shww) # see article provided by Yann
        return min(2000*wind_speed^3,100) # multiply by a scalar for numerical stability. Total cost will be multiplied by a scalar.
        # cap production because for rare data, we get very large values, so we don't really know how to choose the capacity of the substation. 
        # In the future, production scenarios would have to be generated differently.
    end

end

threshold = 1.0 # it is 1.5 officialy, but Yann says it is 1.0 in the reality
scenarios_h = Int.(scenarios_shww .< threshold)
scenarios_p = Production.(scenarios_shww)

# store as text files. One for each strategic period type.
for j in 1:6
    slice_h = scenarios_h[:, j, :]
    writedlm("scenarios_h$(j).txt", slice_h)
    slice_p = scenarios_p[:, j, :]
    writedlm("scenarios_p$(j).txt", slice_p)
end

# code to rebuild scenarios using text files. Used to check we retrieve initial scenarios. 
#scenarios_h_reconstructed = Array{Float64}(undef, 10, 6, 62)
#scenarios_p_reconstructed = Array{Float64}(undef, 10, 6, 62)

#for j in 1:6
#    slice_h = readdlm("scenarios_h$(j).txt")
#    scenarios_h_reconstructed[:, j, :] = slice_h
#    for j in 1:6
#    slice_p = readdlm("scenarios_p$(j).txt")
#    scenarios_p_reconstructed[:, j, :] = slice_p
#    end
#end

