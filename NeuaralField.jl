## Neural Field . Gianni V. Vinci 5/24
include("FP.jl")
using Plots,ProgressMeter,Statistics,DelimitedFiles
theme(:dracula)

# Meta Parameters
NV=150; # Number of grid points
dt=1/10000  # [s];
α=0.25 # Vmin=-αθ
netName="Net.h5"
rateName="rates.dat";

##
cm,S,Npop,Computeμ_σ2=FP.InitializeNet(netName,NV,α,dt)
Life=0.5 # [s]
steps=Int(round(Life/dt))
νD,νN=zeros(Npop),zeros(Npop);
rateFile=open(rateName,"w")


#Simulate
@showprogress for n=1:steps 
    μ,σ2=Computeμ_σ2(νD)
    Threads.@threads for i=1:Npop
        FP.Integrate!(cm[i],μ[i], σ2[i],S[i])
        νD[i]=S[i][cm[i].NV+2]
        νN[i]=S[i][end]
    end
    writedlm(rateFile ,[vcat(n*dt,νN)]) #save output on file
end
close(rateFile);


## Resample and plot average firing rate
t,r=FP.LoadRates(rateName,Life,0.002,dt);
rE=mean(r[:,1:2:end],dims=2)[:,1]
rI=mean(r[:,2:2:end],dims=2)[:,1]
plot(t,rE,label="νₑ")
plot!(t,rI,xlabel="t [s]",ylabel="⟨ν⟩ [Hz]",label="νᵢ")


##savefig("result.pdf")


