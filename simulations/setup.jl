using YAML

p = "input/"
d = YAML.load_file("basecase.yml")

Vs = [0, 1, 2, 3, 4, 5, 10, 50, 100, 1000, -1000, -100, -50, -10, -5, -4, -3, -2, -1.0]
basename = "vseries"

for i = 1:length(Vs)
    outfile = basename * "$i.jld2"
    ymlfile = basename * "$i.yml"
    logfile = basename * "$i.txt"
    d["flow"]["rate"] = 0.45
    if (abs(Vs[i]) > 900) & (abs(Vs[i]) < 4000)
        d["simulation"]["dt"] = 0.00001
    elseif abs(Vs[i]) > 4000
        d["simulation"]["dt"] = 0.000005
    else
        d["simulation"]["dt"] = 0.0001
    end
    d["ions"]["p"] = 1.0e7
    d["potential"]["v1"] = Vs[i]
    d["outputfile"]["name"] = "output/" * outfile
    YAML.write_file(p * ymlfile, d)
end

ions = [16e7, 8e7, 4e7, 2e7, 1e7, 5e6, 2.5e6, 1.25e6]
basename = "ionseries"
p = "input/"
d = YAML.load_file("basecase.yml")
for i = 1:length(ions)
    outfile = basename * "$i.jld2"
    ymlfile = basename * "$i.yml"
    logfile = basename * "$i.txt"
    d["simulation"]["dt"] = 0.00005
    d["ions"]["p"] = ions[i]
    d["potential"]["v1"] = -100
    d["outputfile"]["name"] = "output/" * outfile
    YAML.write_file(p * ymlfile, d)
end

dz = 0:0.1:1.1
basename = "productionseries"
Ls1_0 = d["geometry"]["Ls1"]
Ls2_0 = d["geometry"]["Ls2"]
p = "input/"
d = YAML.load_file("basecase.yml")
for i = 1:length(dz)
    outfile = basename * "$i.jld2"
    ymlfile = basename * "$i.yml"
    logfile = basename * "$i.txt"
    d["geometry"]["Ls1"] = Ls1_0 - dz[i]
    d["geometry"]["Ls2"] = Ls2_0 - dz[i]
    d["simulation"]["dt"] = 0.00005
    d["potential"]["v1"] = -100
    d["outputfile"]["name"] = "output/" * outfile
    YAML.write_file(p * ymlfile, d)
end

Qs = 0.1:0.1:2
basename = "qseries"

for i = 1:length(Qs)
    outfile = basename * "$i.jld2"
    ymlfile = basename * "$i.yml"
    logfile = basename * "$i.txt"
    d["flow"]["rate"] = Qs[i]
    d["ions"]["p"] = 1.0e7
    d["simulation"]["dt"] = 0.0001
    d["potential"]["v1"] = 100.0
    d["outputfile"]["name"] = "output/" * outfile
    YAML.write_file(p * ymlfile, d)
end
