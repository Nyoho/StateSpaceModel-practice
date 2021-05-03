# -*- coding: utf-8 -*-
using StateSpaceModels, CSV, HTTP, DataFrames, Plots
plotly()

# +
using DotEnv
DotEnv.config()

url = ENV["CSV_URL"] # NHKのCOVID-19の都道府県別dailyのCSVファイル
df0 = CSV.File(HTTP.get(url).body) |> DataFrame
# -

function plot_smoothed_state(name)
    df = filter(row -> row[:都道府県名] == name, df0)
    x = convert.(Float64, df.各地の感染者数_1日ごとの発表数)
    t = df.日付
    println("last date: $(t[end])")
    days = 7*8
    xs = x[end-days:end]
    ts = t[end-days:end]
    model = BasicStructural(xs, 7)
    fit!(model)
    print(results(model))
    plt = plot(ts, xs, label="# of cases",
        color="#ccc", shape=:auto,
        marker = (:circle, 4, 0.5, Plots.stroke(1, :gray)),
        xrotation=-90,
        xlabel="date",
        ylabel="number of cases",
        legend=:topleft,
        title="$(name)の最近の$(days)日"
    )
    # filter_output = kalman_filter(model)
    # plot!(plt, ts, get_filtered_state(filter_output), label = "Filtered")
    smoother_output = kalman_smoother(model)
    vs = 2*sqrt.(get_smoothed_state_variance(smoother_output)[1,1,:])
    plot!(plt, ts,
        get_smoothed_state(smoother_output)[:,1],
        ribbon=vs,
        color="dark orange",
        fillalpha=.2,
        ticks=:native,
        label="Smoothed")
end

using Base.Iterators
using Dates
using Printf

function plot_forecast(name; today=-1)
    df = filter(row -> row[:都道府県名] == name, df0)
    x = convert.(Float64, df.各地の感染者数_1日ごとの発表数)
    t = df.日付
    t = map(time_string -> Date(time_string,"y/m/d"), t)
    days = 7*6
    expected_values = []
    σs = []
    for i in reverse(1:days)
        xs = x[end-days+1-i:end-i]
        ts = t[end-days+1-i:end-i]
        model = BasicStructural(xs, 7)
        fit!(model)
        forec = forecast(model, 7)
        e = reduce(vcat,forec.expected_value)
        v = sqrt.(reduce(vcat, forec.covariance))
        push!(expected_values, e[1])
        push!(σs, v[1])
    end
 
    xs = x[end-days+1:end]
    ts = t[end-days+1:end]
    model = BasicStructural(xs, 7)
    fit!(model)
    
    forec = forecast(model, 7)
    latest_e = reduce(vcat,forec.expected_value)
    latest_v = sqrt.(reduce(vcat, forec.covariance))
    
    print("$(name) $(ts[end] + Dates.Day(1))\n")
    @printf("68%%の確率で%.1f〜%.1f人\n", latest_e[1] - latest_v[1], latest_e[1] + latest_v[1])
    @printf("95%%の確率で%.1f〜%.1f人", latest_e[1] - 2*latest_v[1], latest_e[1] + 2*latest_v[1])
    tss = vcat(ts, [ts[end] + Dates.Day(i) for i in 1:7])
    append!(expected_values, latest_e)
    append!(σs, latest_v)
    
    plot(tss, expected_values, ribbon=(2*σs), color="#8ea", fillalpha=.5, label="95% CI")
    #hover_texts = [ [a,b] for (a,b) in zip(expected_values - σs, expected_values + σs) ]
    plot!(tss, expected_values, ribbon=(1*σs), color="#8af", fillalpha=.5, label="68% CI")

    colors = fill("#aaa", length(ts))
    if(today >= 0)
        push!(ts, Dates.today())
        push!(xs, today)
        push!(colors, "#f00")
    end

    plot!(ts, xs, label="# of cases",
        color="#aaa", shape=:auto,
        marker=(:circle, 3, 0.75, Plots.stroke(1, :black)),
        xrotation=-90,
        xlabel="date",
        ylabel="number of cases",
        legend=:topleft,
        ticks=:native,
        title="$(name)の最近の$(days)日と次の7日予測"
    )
    scatter!(ts, xs;
        label="",        
        color=colors,
        marker=(:circle, 3)
    )

end

# 使い方
plot_smoothed_state("広島県")
plot_forecast("東京都")
plot_forecast("東京都", today=759) # 今日のを手動で入れるにはキーワード引数でtodayを設定する。

