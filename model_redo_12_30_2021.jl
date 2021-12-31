@time begin
using StatsBase
using Random
using Distributions
using Plots
using LaTeXStrings

#modification: try to model a process where traits randomly occur in different generations,
# mimicking a type of environmental instability (evolution needs to solve new problems at different generations )
#-----------transfering code from python to here
N=1000
n=7
s1=1
s2=1.5
x1=1
trait_num=18
payoff_spread=0.3

r1=1
r2=1.5
#D=0.6
E=0.0 #payoff bias parameter, depend on the magnitude of s1 and s2 (which is larger)
#when s1>s2, E is positive
w_o_ini=1.0
#variant 1 as the existing variant, starting frequency is 1
ini_T1_perc=0.95
ini_r1_perc=0.95
ini_other1_perc=0.95
beta_ini=0.5
generation=200
mut_mag=0.0

s2_sigma=0.1
r2_sigma=0.1
x2_sigma=0.1

repeat=100

type_mut_rate=0.0

function T1_freq(pop)
    T1_num=0
    for i=1:N
        if pop[i][1]=="T1"
            T1_num=T1_num+1
        end
    end
    return T1_num/N
end

function r1_freq(pop)
    r1_num=0
    for i=1:N
        if pop[i][5]=="r1"
            r1_num=r1_num+1
        end
    end
    return r1_num/N
end

function ave_weight(pop, w)
    w_o_list=[]
    w_a_list=[]
    if w=="w_o"
        for i=1:N
            push!(w_o_list,pop[i][3])
        end
        return(mean(w_o_list))
    elseif w=="w_a"
        for i=1:N
            push!(w_a_list,pop[i][2])
        end
        return(mean(w_a_list))
    end
end

function ave_D(pop)
    D_list=[]
    for i=1:N
        push!(D_list,pop[i][4])
    end
    return mean(D_list)
end

function mutate!(ind)
    ind[2]=ind[2]+rand(Normal(0,mut_mag))
end
#set up initial population
function main(ini_T1_perc,ini_r1_perc,s1,s2,w_o_ini,n,trait_num,r2)

    x2_list=[]
    for i=1:trait_num
        push!(x2_list,rand(Normal(r2, payoff_spread))) #here x2 values are sampled from a normal distribution with mean the same as r2
    end


    function fitness(ind,gen) #include generation in fitness function
        fitness_list=[]
        if ind[1]=="T2"
            push!(fitness_list,rand(Normal(s2, s2_sigma)))
        else
            push!(fitness_list,s1)
        end

        if gen>190 #payoff of the second cultural traits only starts to matter (that is, after they are introduced) after a certain gen
            if ind[5]=="r2" #here, r1/r2 is assigned the last pair of trait
                push!(fitness_list,rand(Normal(r2, r2_sigma)))
            else
                push!(fitness_list,r1)
            end
        end

        for i=6:trait_num+5 #the payoff of these variants only matters after (that is, after they are introduced) a certain generation
            if gen>(i-5)*10
                if ind[i]=="2"
                    push!(fitness_list,rand(Normal(x2_list[i-5],x2_sigma)))
                else
                    push!(fitness_list,x1)
                end
            end
        end
        return (sum(fitness_list))
    end

    pop=[]
    for i=1:N
        if i<=ini_T1_perc*N
            push!(pop,(["T1",rand()*2,w_o_ini,rand()*2-1]))
            #fourth element is D, ranging from -1 to 1
        else
            push!(pop,(["T2",rand()*2,w_o_ini,rand()*2-1]))
        end
    end
    #add more trait (r1 and r2)
    for i=1:N
        if ini_r1_perc>rand()
            push!(pop[i],"r1")
        else
            push!(pop[i],"r2")
        end
    end
    #add even more traits in bulk
    for i=1:N
        for j=1:trait_num
            if ini_other1_perc>rand()
                push!(pop[i],"1")
            else
                push!(pop[i],"2")
            end
        end
    end



    #keep track of T1_freq and population average w_a
    T1_freq_time_series=[]
    r1_freq_time_series=[]

    ave_w_a_time_series=[]
    ave_D_time_series=[]

    for gen=1:generation
        push!(T1_freq_time_series,T1_freq(pop))
        push!(r1_freq_time_series,r1_freq(pop))

        push!(ave_w_a_time_series,ave_weight(pop,"w_a"))
        push!(ave_D_time_series,ave_D(pop))
        fit_list=[]

        for i=1:N
            push!(fit_list,fitness(pop[i],gen))
        end
        #println(fit_list)

        parent_pop=sample(pop,pweights(Vector{Float64}(fit_list)),N,replace=true)
        #here parent_pop are those who have higher fitness and are selected to reproduce.

        #from parent_pop generate a offspring pop, whose w_a mutate
        #here, need to create a offspring_pop from scratch, then copy everything from parent_pop

        offspring_pop=deepcopy(parent_pop)
        for j=1:N
            mutate!(offspring_pop[j])#mutate weight, change of offspring_pop doesn't affect parent_pop
        end


        #now each offspring sample n number of models from pop at random

        for k=1:N #k is index of offspring

            model_sample=sample(pop,n)

            #compute the number of T1 indivdiuals in model set:
            T1_num=0
            T2_num=0
            for model_indx=1:n
                if model_sample[model_indx][1]=="T1"
                    T1_num=T1_num+1
                elseif model_sample[model_indx][1]=="T2"
                    T2_num=T2_num+1
                end
            end

            #compute the number of T1 indivdiuals in model set:
            T1_num=0
            T2_num=0
            for model_indx=1:n
                if model_sample[model_indx][1]=="T1"
                    T1_num=T1_num+1
                elseif model_sample[model_indx][1]=="T2"
                    T2_num=T2_num+1
                end
            end


            if T1_num==0
                p_T1_adopted=0

            elseif T2_num==0
                p_T1_adopted=1

            elseif T1_num>n/2
                p_T1_adopted=((T1_num+offspring_pop[k][4])*offspring_pop[k][2]+(s1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(s1+s2)*offspring_pop[k][3])

            elseif T1_num==n/2
                p_T1_adopted=((T1_num)*offspring_pop[k][2]+(s1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(s1+s2)*offspring_pop[k][3])

            elseif T1_num<n/2
                p_T1_adopted=((T1_num-offspring_pop[k][4])*offspring_pop[k][2]+(s1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(s1+s2)*offspring_pop[k][3])

            end
            #println(p_T1_adopted)
            if p_T1_adopted>rand()
                offspring_pop[k][1]="T1"
            else
                offspring_pop[k][1]="T2"
            end

            #now do the same thing for r1/r2:
            #compute the number of r1 indivdiuals in model set:
            r1_num=0
            r2_num=0
            for model_indx=1:n
                if model_sample[model_indx][5]=="r1"
                    r1_num=r1_num+1
                elseif model_sample[model_indx][5]=="r2"
                    r2_num=r2_num+1
                end
            end
            if r1_num==0
                p_r1_adopted=0

            elseif r2_num==0
                p_r1_adopted=1

            elseif r1_num>n/2
                p_r1_adopted=((r1_num+offspring_pop[k][4])*offspring_pop[k][2]+(r1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(r1+r2)*offspring_pop[k][3])

            elseif r1_num==n/2
                p_r1_adopted=((r1_num)*offspring_pop[k][2]+(r1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(r1+r2)*offspring_pop[k][3])

            elseif r1_num<n/2
                p_r1_adopted=((r1_num-offspring_pop[k][4])*offspring_pop[k][2]+(r1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(r1+r2)*offspring_pop[k][3])

            end
            if gen>190 #p_r1_adopted only matters when gen>100
                #println(p_T1_adopted)
                if p_r1_adopted>rand()
                    offspring_pop[k][5]="r1"
                else
                    offspring_pop[k][5]="r2"
                end
            else #before 10 generations
                if ini_r1_perc>rand()
                    offspring_pop[k][5]="r1"
                else
                    offspring_pop[k][5]="r2"
                end
            end

            #now do the same thing for the rest of traits
            for x=6:trait_num+5

                x1_num=0
                x2_num=0
                for model_indx=1:n
                    if model_sample[model_indx][x]=="1"
                        x1_num=x1_num+1
                    elseif model_sample[model_indx][x]=="2"
                        x2_num=x2_num+1
                    end
                end
                if x1_num==0
                    p_x1_adopted=0

                elseif x2_num==0
                    p_x1_adopted=1

                elseif x1_num>n/2
                    p_x1_adopted=((x1_num+offspring_pop[k][4])*offspring_pop[k][2]+(x1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(x1+x2_list[x-5])*offspring_pop[k][3])

                elseif x1_num==n/2
                    p_x1_adopted=((x1_num)*offspring_pop[k][2]+(x1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(x1+x2_list[x-5])*offspring_pop[k][3])

                elseif x1_num<n/2
                    p_x1_adopted=((x1_num-offspring_pop[k][4])*offspring_pop[k][2]+(x1+E)*offspring_pop[k][3])/(n*offspring_pop[k][2]+(x1+x2_list[x-5])*offspring_pop[k][3])
                end

                if gen>(x-5)*10 #p_r1_adopted only matters when gen is larger than a certain number
                    if p_r1_adopted>rand()
                        offspring_pop[k][x]="1"
                    else
                        offspring_pop[k][x]="2"
                    end
                else
                    if ini_other1_perc>rand()
                        offspring_pop[k][x]="1"
                    else
                        offspring_pop[k][x]="2"
                    end
                end
            end

        end

    pop=offspring_pop
    end

    global T1_freq_time_series
    global r1_freq_time_series
    global ave_w_a_time_series
    global ave_D_time_series
    p1=plot((1:generation),T1_freq_time_series, ylabel="T1 frequency")
    ylims!((0.0,1.0))
    plot!([s1/(s1+s2)],seriestype="hline")

    p1_5=plot((1:generation),r1_freq_time_series, ylabel="r1 frequency")
    ylims!((0.0,1.0))
    plot!([r1/(r1+r2)],seriestype="hline")


    p2=plot((1:generation),ave_w_a_time_series,xlabel="generation", ylabel=L"w_a")
    ylims!((0.0,2.0))

    p3=plot((1:generation),ave_D_time_series,xlabel="generation", ylabel=L"D")
    ylims!((-0.5,1.0))
    plot(p1, p1_5, p2, p3,layout = (4, 1), legend = false)
end




#create same timeseries over multiple runs
s2_r2_set=[0.8,1.2]
weight_array=[]
D_array=[]
freq_array=[]
freq_array_r=[]

for i in s2_r2_set
    repeat_set_weight=[]
    repeat_set_D=[]
    repeat_set_freq=[]
    repeat_set_freq_r=[]
    for j in 1:repeat
        main(ini_T1_perc,ini_r1_perc,s1,i,w_o_ini,n,trait_num,i)
        push!(repeat_set_freq, T1_freq_time_series)
        push!(repeat_set_freq_r,r1_freq_time_series)
        push!(repeat_set_weight, ave_w_a_time_series)
        push!(repeat_set_D,ave_D_time_series)
        println("s2_r2= ", i, "repeat= ", j)
    end
    weight_ave=mean(repeat_set_weight, dims=1)
    D_ave=mean(repeat_set_D,dims=1)
    freq_ave=mean(repeat_set_freq,dims=1)
    freq_ave_r=mean(repeat_set_freq_r,dims=1)
    push!(weight_array, weight_ave)
    push!(D_array,D_ave)
    push!(freq_array, freq_ave)
    push!(freq_array_r, freq_ave_r)
end

p_D_combined=plot((1:generation),freq_array,label = [L"s_2=x_2=0.6" L"s_2=x_2=1.4" ], ylabel="C1 freq",
xtickfontsize=11,ytickfontsize=10,xguidefontsize=13,yguidefontsize=12,color=[:black :purple :red])
ylims!((-0.1,1.1))
#plot!([0.8/(1+0.8)],seriestype="hline",label="y=s1/(s1+s2)",line = (:dot,2),color=[:black])
#plot!([1.2/(1+1.2)],seriestype="hline",label="y=r1/(r1+r2)",line = (:dot,2),color=[:black])
r_D_combined=plot((1:generation),freq_array_r,label = ["trait_num=5" "trait_num=10"], ylabel="R1 freq",
xtickfontsize=11,ytickfontsize=10,xguidefontsize=13,yguidefontsize=12,color=[:black :purple :red])
ylims!((-0.1,1.1))
#plot!([1.2/(1+1.2)],seriestype="hline",label="y=r1/(r1+r2)",line = (:dot,2),color=[:black])

w_D_combined=plot((1:generation),weight_array,label=["trait_num=5" "trait_num=10"], ylabel=L"w_a",
xtickfontsize=11,ytickfontsize=10,xguidefontsize=13,yguidefontsize=20,color=[:black :purple :red])
ylims!((0.0,2))
plot!([1],seriestype="hline",label="reference line",line = (:dot,2),color=[:black])

D_D_combined=plot((1:generation),D_array,label=["trait_num=5" "trait_num=10"],xlabel="generation", ylabel=L"D",
xtickfontsize=11,ytickfontsize=10,xguidefontsize=13,yguidefontsize=20,color=[:black :purple :red])
ylims!((-0.2,0.2))
plot!([0],seriestype="hline",label="reference line D=0",line = (:dot,2),color=[:black])

plot(p_D_combined, w_D_combined, D_D_combined,layout = (3,1), legend = false,size = (250, 500))

savefig("8_4_new.png") #saved to Documents folder
end#final, for time
