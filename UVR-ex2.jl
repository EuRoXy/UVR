### A Pluto.jl notebook ###
# v0.12.15

using Markdown
using InteractiveUtils

# ╔═╡ d2edd716-318d-11eb-03c7-5d60affb3336
using Statistics

# ╔═╡ b4d32cdc-318c-11eb-2a4b-290dbd516c5b
using DataFrames

# ╔═╡ 1b7ffdb8-3258-11eb-0617-e5e1dd0b4610
using DelimitedFiles

# ╔═╡ 6461ecf4-3261-11eb-3f8f-8f287d22188d
using NumericalIntegration

# ╔═╡ 1287f4f2-3494-11eb-38b9-a3ba81f6d0b4
using Plots

# ╔═╡ 331c6c28-3176-11eb-092c-adba2db0abd7
md"
### Exercise 2 UV Radiation
##### 1)
Assuming we have ozone hole conditions over middle latitudes, calculate the **UV Index** in Zurich at local **noon on 1 June** 2020 

assuming clear sky conditions for total column **ozone of 94 DU** (minimum of 2020 antarctic ozone hole) and 

compare the results to standard conditions at Zurich.
"

# ╔═╡ 0f0e87f2-318b-11eb-3070-1db2337e96b8
# solar elevation angle by https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
δ = 59.24

# ╔═╡ a73c1e68-318b-11eb-2df9-793449b23bab
# solar zenith angle
θ = 90 - δ

# ╔═╡ cc606208-318b-11eb-3565-53093ce68aa0
md"in1.inp
```
data_files_path /opt/libRadtran/data/
atmosphere_file /opt/libRadtran/data/atmmod/afglus.dat
source solar /opt/libRadtran/data/solar_flux/atlas_plus_modtran
mol_modify O3 94. DU # Set ozone column
day_of_year 153 # Correct for Earth- Sun distance
albedo 0.2 # Surface albedo
sza 30.8 # Solar zenith angle
rte_solver disort # Radiative transfer equation solver
number_of_streams 6 # Number of streams
wavelength 299.0 341.0 # Wavelength range [nm]
```"

# ╔═╡ 1e3119ae-3190-11eb-2778-b9a8919d680b
md"Assume standard total ozone column over Zurich to be 310 DU.

ref: [MeteoSwiss](https://www.meteoswiss.admin.ch/home/climate/climate-change-in-switzerland/ozone-monitoring.html)"

# ╔═╡ d90d536e-317f-11eb-2cb9-dd1de80b495f
md"in2.inp
```
mol_modify O3 310 DU # Set ozone column
```"

# ╔═╡ d0db06e6-3193-11eb-0d86-a7c896d1af75
sp1 = DataFrame(readdlm("out1.csv", '\t'), ["Wavelength [nm]", "Direct beam irradiance [mW m⁻² nm⁻¹]", "Diffuse down irradiance [mW m⁻² nm⁻¹]"]);

# ╔═╡ ddf93e70-325f-11eb-357f-3ff1bc7f8e9b
begin
	sp1["Total irradiance [W m⁻² nm⁻¹]"] = (sp1["Direct beam irradiance [mW m⁻² nm⁻¹]"] .+ sp1["Diffuse down irradiance [mW m⁻² nm⁻¹]"]) ./ 1000
	sp1
end

# ╔═╡ 0d5b89ce-3261-11eb-2ce7-7f1294efadc9
begin
	sp2 = DataFrame(readdlm("out2.csv", '\t'), ["Wavelength [nm]", "Direct beam irradiance [mW m⁻² nm⁻¹]", "Diffuse down irradiance [mW m⁻² nm⁻¹]"])
	sp2["Total irradiance [W m⁻² nm⁻¹]"] = (sp2["Direct beam irradiance [mW m⁻² nm⁻¹]"] .+ sp2["Diffuse down irradiance [mW m⁻² nm⁻¹]"]) ./ 1000
end

# ╔═╡ d29bc9d4-3260-11eb-2dca-4fd1fd702454
function heaviside(t)
   0.5 * (sign(t) + 1)
end

# ╔═╡ db305178-3260-11eb-27a2-2386fc933704
function interval(t, a, b)
   heaviside(t-a) - heaviside(t-b)
end

# ╔═╡ e2306c40-3260-11eb-170e-f55ef8653cba
function erythemaAS(λ)
    1. .* interval(λ, 250,298) + 
    10^(0.094*(298-λ)) .* interval(λ, 298,328) + 
    10^(0.015*(140-λ)) .* interval(λ, 328,400)
end

# ╔═╡ 27c364a8-3261-11eb-245e-95ff094dfa65
begin
	eryWeiIrr1 = sp1["Total irradiance [W m⁻² nm⁻¹]"] .* erythemaAS.(sp1["Wavelength [nm]"])
	Eery1 = integrate(sp1["Wavelength [nm]"], eryWeiIrr1)
end

# ╔═╡ 712828e0-3261-11eb-0766-f960029f05eb
begin
	eryWeiIrr2 = sp2["Total irradiance [W m⁻² nm⁻¹]"] .* erythemaAS.(sp2["Wavelength [nm]"])
	Eery2 = integrate(sp2["Wavelength [nm]"], eryWeiIrr2)
end

# ╔═╡ 8f208324-3261-11eb-0a1e-0bfed9111991
UVI1, UVI2 = round(Eery1*40, digits=2), round(Eery2*40, digits=2)

# ╔═╡ a840f944-3261-11eb-33e6-591a610f8522
md"_The UV index in the first situation is **$UVI1**, which is around 3 times of the UV index of **$UVI2** in the standard situation._"

# ╔═╡ 2ef9cecc-3489-11eb-265d-156bf386b47a
md"
##### 2)
Estimate the impact of the California Forest fires on the erythemal weighted irradiance on 10 September. 
Use aerosol info from [Aeronet site Monterey](https://aeronet.gsfc.nasa.gov/cgi-bin/data_display_aod_v3?site=Monterey&nachal=0&year=2020&month=9&day=10&aero_water=0&level=1&if_day=0&if_err=0&place_code=10&year_or_month=0)
"

# ╔═╡ bd2ffaf8-3499-11eb-146e-d1d524ba1eb2
md"in30.inp
```
data_files_path /opt/libRadtran/data/
atmosphere_file /opt/libRadtran/data/atmmod/afglus.dat
source solar /opt/libRadtran/data/solar_flux/atlas_plus_modtran

day_of_year 254 # Correct for Earth- Sun distance
albedo 0.2 # Surface albedo
sza 35.3 # Solar zenith angle at 19:00 UTC
rte_solver disort # Radiative transfer equation solver
number_of_streams 6 # Number of streams
wavelength 299.0 341.0 # Wavelength range [nm]
```"

# ╔═╡ a6a6f03a-3489-11eb-227c-678174de740e
md"in3.inp
```
aerosol_vulcan 1 # Aerosol type above 2km
aerosol_haze 6 # Aerosol type below 2km
aerosol_season 2 # Fall-winter profile.
aerosol_visibility 20.0 # Visibility
aerosol_modify tau set 2.8
```"

# ╔═╡ 04337414-349c-11eb-1e7b-8dcd0fac751b
md"_California forest fires reduced the the erythemal weighted irradiance on 10 September to **55% of its normal level**._"

# ╔═╡ 5978e5d4-3489-11eb-25ba-f9cd1868336f
md"
##### 3)
Calculate the spectral irradiance (290 to 400 nm) for several albedo values between 0.04 (grass), and 0.9 (fresh snow) and 

two total column ozone amounts (for example 300 and 500 DU). 

Plot graphically (ratio relative to grass albedo) and discuss the results.
"

# ╔═╡ d5cb888e-348a-11eb-3c24-f7bd592447db
md"in4.inp
```
data_files_path /opt/libRadtran/data/
atmosphere_file /opt/libRadtran/data/atmmod/afglus.dat
source solar /opt/libRadtran/data/solar_flux/atlas_plus_modtran
mol_modify O3 300. DU # Set ozone column
day_of_year 170 # Correct for Earth- Sun distance
albedo 0.04 # surface albedo
sza 32 # Solar zenith angle
rte_solver disort # Radiative transfer equation solver
number_of_streams 6 # Number of streams
wavelength 290.0 400.0 # Wavelength range [nm]
```"

# ╔═╡ f1acd45e-348a-11eb-1b5c-35710d110ca7
md"in5.inp
```
albedo 0.4 # Surface albedo
```"

# ╔═╡ f7446b80-348e-11eb-02de-d53b7b78067c
md"in6.inp
```
albedo 0.9 # Surface albedo
```
"

# ╔═╡ 68062e58-348f-11eb-00f4-b757c2e4f404
md"in7.inp
```
mol_modify O3 400. DU # Set ozone column
albedo 0.04 # surface albedo
```"

# ╔═╡ 84f1a560-348f-11eb-19b4-81d9806da108
md"in8.inp
```
mol_modify O3 500. DU # Set ozone column
albedo 0.04 # surface albedo
```"

# ╔═╡ 669cb44e-3497-11eb-1db4-0fa34f813f35
function spec(fnam)
	sp = DataFrame(readdlm(fnam, '\t'), ["Wavelength [nm]", "Direct beam irradiance [mW m⁻² nm⁻¹]", "Diffuse down irradiance [mW m⁻² nm⁻¹]"])
	sp["Total irradiance [W m⁻² nm⁻¹]"] = (sp["Direct beam irradiance [mW m⁻² nm⁻¹]"] .+ sp["Diffuse down irradiance [mW m⁻² nm⁻¹]"]) ./ 1000
	x = sp["Wavelength [nm]"]
	y = sp["Total irradiance [W m⁻² nm⁻¹]"]
	return x, y
end

# ╔═╡ 937e6c0c-349a-11eb-1e71-01d3911cb92e
begin
	sp30 = spec("out30.csv");
	eryWeiIrr30 = sp30[2] .* erythemaAS.(sp30[1])
	Eery30 = integrate(sp30[1], eryWeiIrr30)
end

# ╔═╡ e02733e0-349a-11eb-0936-ad198ad05739
begin
	sp3 = spec("out3.csv");
	eryWeiIrr3 = sp3[2] .* erythemaAS.(sp3[1])
	Eery3 = integrate(sp3[1], eryWeiIrr3)
end

# ╔═╡ f19510fc-349a-11eb-12d2-dbcd59b0bc6c
imp = Eery3 / Eery30

# ╔═╡ d46376d6-3497-11eb-25e8-a7183afeb207
begin
	sp4 = spec("out4.csv");
	sp5 = spec("out5.csv");
	sp6 = spec("out6.csv");
	sp7 = spec("out7.csv");
	sp8 = spec("out8.csv");
end

# ╔═╡ a9c3da6a-3497-11eb-092c-930f36809a99
begin
	p1 = plot(sp4[1], sp4[2], label="α=0.04", leg=:topleft, 
			  title="O3 = 300DU", labelfontsize=8,
		      xlabel="Wavelength [nm]", ylabel="Irradiance [W m⁻² nm⁻¹]", )
	plot!(sp5[1], sp5[2], label="α=0.4")
	plot!(sp6[1], sp6[2], label="α=0.9")
end

# ╔═╡ 809969b8-349d-11eb-273a-9d54e9675bb7
begin
	p3 = plot(sp4[1], sp4[2]./sp4[2], label="α=0.04", leg=:topright, xlim=(290,400),
			  title="O3 = 300DU", labelfontsize=8, xlabel="Wavelength [nm]", 					  ylabel="Ralative Ratio of spectral irradiance", )
	plot!(sp5[2]./sp4[2], label="α=0.4")
	plot!(sp6[2]./sp4[2], label="α=0.9")
end

# ╔═╡ f6ea15f0-349c-11eb-1a65-afce8aafeef0
md"The spectral irradiance generally **increases with increasing albedo** given the same ozone concentration, the wavelength dependence is small;

The spectral irradiance also **decreases with increasing ozone concentration** given the same albedo, the relative ratio increases with increasing wavelength.
"

# ╔═╡ e3393dc8-3498-11eb-3e60-a1ccf7a4aa88
begin
	p2 = plot(sp4[1], sp4[2], label="O3=300DU", leg=:topleft, title="α = 0.04",
		 xlabel="Wavelength [nm]", ylabel="Irradiance [W m⁻² nm⁻¹]", labelfontsize=8)
	plot!(sp7[1], sp7[2], label="O3=400DU")
	plot!(sp8[1], sp8[2], label="O3=500DU")

end

# ╔═╡ 2ef63a04-349e-11eb-309f-dd5b98483a02
begin
	p4 = plot(sp4[1], sp4[2]./sp4[2], label="O3=300DU", leg=:topleft, 
			  title="α = 0.04", labelfontsize=8, xlim=(290,400),
			  xlabel="Wavelength [nm]",ylabel="Ralative Ratio of spectral irradiance")
	plot!(sp7[2]./sp4[2], label="O3=400DU")
	plot!(sp8[2]./sp4[2], label="O3=500DU")
end

# ╔═╡ Cell order:
# ╟─331c6c28-3176-11eb-092c-adba2db0abd7
# ╠═0f0e87f2-318b-11eb-3070-1db2337e96b8
# ╠═a73c1e68-318b-11eb-2df9-793449b23bab
# ╠═d2edd716-318d-11eb-03c7-5d60affb3336
# ╠═b4d32cdc-318c-11eb-2a4b-290dbd516c5b
# ╠═1b7ffdb8-3258-11eb-0617-e5e1dd0b4610
# ╠═6461ecf4-3261-11eb-3f8f-8f287d22188d
# ╠═1287f4f2-3494-11eb-38b9-a3ba81f6d0b4
# ╟─cc606208-318b-11eb-3565-53093ce68aa0
# ╟─1e3119ae-3190-11eb-2778-b9a8919d680b
# ╟─d90d536e-317f-11eb-2cb9-dd1de80b495f
# ╠═d0db06e6-3193-11eb-0d86-a7c896d1af75
# ╠═ddf93e70-325f-11eb-357f-3ff1bc7f8e9b
# ╠═0d5b89ce-3261-11eb-2ce7-7f1294efadc9
# ╠═d29bc9d4-3260-11eb-2dca-4fd1fd702454
# ╠═db305178-3260-11eb-27a2-2386fc933704
# ╠═e2306c40-3260-11eb-170e-f55ef8653cba
# ╠═27c364a8-3261-11eb-245e-95ff094dfa65
# ╠═712828e0-3261-11eb-0766-f960029f05eb
# ╠═8f208324-3261-11eb-0a1e-0bfed9111991
# ╟─a840f944-3261-11eb-33e6-591a610f8522
# ╟─2ef9cecc-3489-11eb-265d-156bf386b47a
# ╟─bd2ffaf8-3499-11eb-146e-d1d524ba1eb2
# ╟─a6a6f03a-3489-11eb-227c-678174de740e
# ╠═937e6c0c-349a-11eb-1e71-01d3911cb92e
# ╠═e02733e0-349a-11eb-0936-ad198ad05739
# ╠═f19510fc-349a-11eb-12d2-dbcd59b0bc6c
# ╟─04337414-349c-11eb-1e7b-8dcd0fac751b
# ╟─5978e5d4-3489-11eb-25ba-f9cd1868336f
# ╟─d5cb888e-348a-11eb-3c24-f7bd592447db
# ╟─f1acd45e-348a-11eb-1b5c-35710d110ca7
# ╟─f7446b80-348e-11eb-02de-d53b7b78067c
# ╟─68062e58-348f-11eb-00f4-b757c2e4f404
# ╟─84f1a560-348f-11eb-19b4-81d9806da108
# ╠═669cb44e-3497-11eb-1db4-0fa34f813f35
# ╠═d46376d6-3497-11eb-25e8-a7183afeb207
# ╠═a9c3da6a-3497-11eb-092c-930f36809a99
# ╠═809969b8-349d-11eb-273a-9d54e9675bb7
# ╟─f6ea15f0-349c-11eb-1a65-afce8aafeef0
# ╠═e3393dc8-3498-11eb-3e60-a1ccf7a4aa88
# ╠═2ef63a04-349e-11eb-309f-dd5b98483a02
