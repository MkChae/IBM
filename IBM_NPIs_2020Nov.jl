""" 
    2023-08-01
    IBM Julia language code developed by Min-Kyung Chae (in NIMS).
    Multi-thread code:     julia --threads 2 filename.jl

    Julia version: 19

    Individual-based model (IBM) to simulate the spread of COVID-19 in South Korea

    Intervention policy : November 2020
    Level 1 : school capacity reduce + large workplace telework
    Level 2 : private gathering size is maximum 4 in the national capital area
    Level 3 : private gathering size is maximum 4 in the entire country
"""

using Base.Threads 
using DataFrames, CSV, StatsBase, Statistics, Distributions # DataFrames(v1.6.1), CSV(v0.10.11), StatsBase(v0.33.21), Statistics(v1.9.0), Distributions(v0.25.98)
using Random
using Plots, Images # Plots(v1.38.17), Images(v0.25.3)
using RData # v1.0.0
using NaNMath # v1.0.2
import CodecBzip2 # v0.7.2


function Relative_Infectiousness(std_, rngs_)

    if std_ == 0
        rn = 1
    else
        mean_ = 1;
        k = (mean_*mean_)/(std_*std_);     theta = (std_*std_)/mean_;
        rn = round(rand(rngs_, Gamma(k, theta)), digits=3)
    end

    return rn
end

function InitialSetting_vs()
    vs_ = zeros(10)
    vs_[1]  = 0.12956;    vs_[2]  = 0.236426;    vs_[3]  = 0.213072;    vs_[4]  = 0.156682;    vs_[5]  = 0.104453;
    vs_[6]  = 0.0657157;  vs_[7]  = 0.0397921;   vs_[8]  = 0.0234487;   vs_[9]  = 0.0135399;   vs_[10] = 0.00769614;

    return vs_
end

function State_E_Period(rngs_)
    k = 1.926;     theta = 1.775;
    rn = round(Int64, rand(rngs_, Gamma(k, theta)))
    rn = (rn < 1 ? 1 : rn)
    rn = (rn > 10 ? 10 : rn)
    
    return rn
end

function State_I_Period()
    rn = 8
    
    return rn  
end

function Friends_Meeting_Number(n_, rngs_)
    out_max = 9  # Meeting size should be +1. including me. (i.e. max size 10)
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Friends_Meeting_Number_Constraint(df_C_, rc_, n_, rngs_)
    out_max = df_C_[rc_, :constraint] - 1  # Meeting size should be +1. including me.
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Friends_Meeting_Number_Constraint_Lv3(n_, rngs_)
    out_max = 3  # Meeting size should be +1. including me. (i.e. max size 4)
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Setting_Group_Education_Rotation(df_, group1_, group2_, group3_, group4_, rngs_)
    # Only 30% of the classes in each region go to school (1/3 of them go to school) > They are assigned a group to go to school.
    # It's very time consuming because we have to go through 250 region, but we only need to do it once.

    k = 1;     school_group_ = [];
    for cols in 4:7 # Kindergarten, elementary, junior, and high school
        for i in 1:250
            list = unique(df_[findall((df_[:, :region_code] .== i) .& (df_[:, cols] .!= 0) .& (df_[:, :age] .<= 18)), cols])
            list_rn = zeros(length(list))
            list_rn .= [rand(rngs_) for i in 1:length(list)]

            if i == 1
                for j in 1:3
                    if cols == 7
                        push!(school_group_, symdiff(list, list[findall((j-1.0)/3.0 .<= list_rn .< j/3.0)]))
                    else
                        push!(school_group_, list[findall((j-1.0)/3.0 .<= list_rn .< j/3.0)])
                    end
                end
            else
                if cols == 7
                    append!(school_group_[k], symdiff(list, list[findall(0.0/3.0 .<= list_rn .< 1.0/3.0)]))
                    append!(school_group_[k+1], symdiff(list, list[findall(1.0/3.0 .<= list_rn .< 2.0/3.0)]))
                    append!(school_group_[k+2], symdiff(list, list[findall(2.0/3.0 .<= list_rn .< 3.0/3.0)]))
                else
                    append!(school_group_[k], list[findall(0.0/3.0 .<= list_rn .< 1.0/3.0)])
                    append!(school_group_[k+1], list[findall(1.0/3.0 .<= list_rn .< 2.0/3.0)])
                    append!(school_group_[k+2], list[findall(2.0/3.0 .<= list_rn .< 3.0/3.0)])
                end
            end
        end
        
        k += 3
    end

    k = 1;    tmp0 = [length(group1_), length(group2_), length(group3_), length(group4_)]
    for i in 1:4
        tmp = length(unique(union(school_group_[k], school_group_[k+1], school_group_[k+2])))
        if tmp+1 != tmp0[i]
            println(" ??? ERROR in setting group education rotation ", tmp, " ", tmp0[i])
        end
        k += 3
    end

    group1_ = nothing; group2_ = nothing; group3_ = nothing; group4_ = nothing;
    tmp0 = nothing; list = nothing; list2 = nothing;
    GC.gc()

    return school_group_
end

function ReadCSVFile()
    df1 = DataFrame(CSV.File("./src/synthetic_population_info.csv"))
    df2 = DataFrame(CSV.File("./src/synthetic_population_friends.csv"))
    df3 = DataFrame(CSV.File("./src/synthetic_population_religion.csv"))
    df4 = DataFrame(CSV.File("./src/synthetic_population_agegroup.csv"))
    df5 = DataFrame(CSV.File("./src/synthetic_population_groupnumber.csv"))
    
    df_ = outerjoin(df1, df2; on=:person_id)
    df_ = outerjoin(df_, df3; on=:person_id)
    df_ = outerjoin(df_, df4; on=:person_id)
    df_ = outerjoin(df_, df5; on=:person_id)

    df1 = nothing; df2 = nothing; df3 = nothing; df4 = nothing; df5 = nothing;
    GC.gc()

    return df_
end

function InfodfSetting(df_)
    df_[!, 1:10] = convert.(Int64, df_[:, 1:10])
    df_[!, 12:19] = convert.(Int64, df_[:, 12:19])

    df_[!, :"tmp"] .= [Vector{Int64}()]
    for i in 1:nrow(df_)
        if df_[i, :n_friends] != 0
            list = parse.(Int, split(chop(df_[i, :friends]; head=1, tail=1), ','))
            df_[i, :tmp] = list
            list = nothing
        end
    end
    select!(df_, Not(:friends))
    rename!(df_, :tmp => :friends)

    return df_
end

function MakedfDay(N_)
    df_D_ = DataFrame("person_id" => 1:N_, 
    "religion_out" => 0,  
    "friends_out" => 0, "friends_m_num" => 0, 
    "now_region" => 0, "encounter_out" => 0, "encounter_m_num" => 0)

    return df_D_
end

function MakedfState(N_)
    df_S_ = DataFrame("person_id" => 1:N_, "state" => 0,  
    "period_E" => 0, "period_I" => 0, "date_E" => 0, "date_I" => 0, "date_R" => 0, 
    "date_unt" => 0, "date_inf" => 0, "vs_count" => 0, 
    "day_count" => 0)
    df_S_[!, :"I_type"] .= [Vector{Int64}()]
    df_S_[!, :"I_pop"] .= [Vector{Int64}()]
    df_S_[!, :"RI"] .= 0.0
    df_S_[!, :"vs_N"] .= 0.0

    return df_S_
end

function MakedfLambda(N_)
    df_I_ = DataFrame("person_id" => 1:N_)
    df_I_[!, :"I_pop_h"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_h"] .= 0
    df_I_[!, :"N_h"] .= 0
    df_I_[!, :"I_pop_s"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_s"] .= 0
    df_I_[!, :"N_s"] .= 0
    df_I_[!, :"I_pop_o"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_o"] .= 0
    df_I_[!, :"N_o"] .= 0
    df_I_[!, :"I_pop_r"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_r"] .= 0
    df_I_[!, :"N_r"] .= 0
    df_I_[!, :"I_pop_f"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_f"] .= 0
    df_I_[!, :"N_f"] .= 0
    df_I_[!, :"I_pop_e"] .= [Vector{Int64}()]
    df_I_[!, :"N_pop_e"] .= 0
    df_I_[!, :"N_e"] .= 0

    return df_I_
end

# Meeting restrictions data frame
function MakeRegionConstraint()
    df_C_ = DataFrame("region_code" => 1:250, "constraint" => 10, "beta" => 1.0) 
    # Up to 4 individuals in the national capital area (10 elsewhere)
    df_C_[1:25, :constraint] .= 4; df_C_[1:25, :beta] .= 0.5;  # Seoul
    df_C_[50:59, :constraint] .= 4; df_C_[50:59, :beta] .= 0.5;  # Incheon
    df_C_[75:116, :constraint] .= 4; df_C_[75:116, :beta] .= 0.5; # Gyeonggi

    return df_C_
end

function MovePopulation(group_region_, rngs_)
    other_region_list_ = []

    TR_OD = DataFrame(CSV.File("./src/train.csv"))
    PL_OD = DataFrame(CSV.File("./src/plain.csv"))
    
    for rc in 1:250
        tmp = group_region_[rc]
        n0 = nrow(tmp)
        df_tmp = DataFrame("id" => 1:n0, "check" => 0)
        list = zeros(Int64, n0)
        
        for j in 1:250
            # train
            id = df_tmp[(df_tmp.check .== 0), :id]
            n = length(id)
            sample_rows = sample(rngs_, 1:n, round(Int64, n*TR_OD[rc, j]), replace=false)
            df_tmp[id[sample_rows], :check] .= 1
            list[id[sample_rows]] .= j
    
            # plane
            id = df_tmp[(df_tmp.check .== 0), :id]
            n = length(id)
            sample_rows = sample(rngs_, 1:n, round(Int64, n*PL_OD[rc, j]), replace=false)
            df_tmp[id[sample_rows], :check] .= 1
            list[id[sample_rows]] .= j
        end
    
        id = df_tmp[(df_tmp.check .== 0), :id]
        list[id] .= rc
        push!(other_region_list_, list)

    end
    
    TR_OD = nothing; PL_OD = nothing;
    tmp = nothing; df_tmp = nothing;
    id = nothing; list = nothing; sample_rows = nothing; id = nothing;
    group_region_ = nothing;
    GC.gc()

    return other_region_list_
end

# Print person ID lists who in the household has been in contact with the infected person
function House_Check(I_list_, df_, df_S_, df_I_, group_house_)

    # 1 Output the house_id of the household with the infected person
    contact_house_id = unique(df_[I_list_, :house_id])
    contact_house_N = length(contact_house_id)
    
    # 2 Print the person_id of each household member living in the household with the infected person.
    contact_house_person_id = [group_house_[contact_house_id[i]][:, :person_id] for i in 1:contact_house_N]
    
    # 3 Make person_id a one-dimensional vector (to make it easier to add contact_list later)
    contact_house_person_id_1D = Vector{Int64}()
    [append!(contact_house_person_id_1D, contact_house_person_id[i]) for i in 1:contact_house_N]
    
    # 4 Print the status of each household member
    contact_house_person_state = [df_S_[contact_house_person_id[i], :state] for i in 1:contact_house_N]
    
    # 5 Record the number of infected people in the household where each member lives.
    contact_house_I_N = [length(contact_house_person_state[i][contact_house_person_state[i] .== 2]) for i in 1:contact_house_N]
    
    # 6 Record information about each person in df_I for lambda computation
    # N_pop_h : Number of people infected in the household, I_pop_h : person id of the infected person in the household, N_h : Household size
    for i in 1:contact_house_N
        df_I_[contact_house_person_id[i], :N_pop_h] .= contact_house_I_N[i]
        df_I_[contact_house_person_id[i], :I_pop_h] .= [contact_house_person_id[i][findall(x -> x == 2, df_S_[contact_house_person_id[i], :state])] for j in 1:length(contact_house_person_id[i])]
        df_I_[contact_house_person_id[i], :N_h] .= length(contact_house_person_id[i])
    end
    
    list = contact_house_person_id_1D[findall(x -> x == 0, df_S_[contact_house_person_id_1D, :state])]

    group_house_ = nothing; I_list_ = nothing;
    contact_house_id = nothing; contact_house_person_id = nothing; contact_house_person_id_1D = nothing;
    contact_house_person_state = nothing; contact_house_I_N = nothing;

    return list
end

# Print out who in your office has been in contact with an infected person
function Office_Check(I_list_, df_, df_S_, df_I_, group_office_)

    # 0 Print out the employees who are infected.
    I_list_office = I_list_[findall(x -> x != 0, df_[I_list_, :office])]
    
    if length(I_list_office) != 0
        # 1 Output the office_num of the office where the infected person is located.
        contact_office_id = unique(df_[I_list_office, :office])
        contact_office_N = length(contact_office_id)

        # 2 Output the person_id of the members who work in the office where the infected person is located.
        contact_office_person_id = [group_office_[contact_office_id[i]+1][:, :person_id] for i in 1:contact_office_N]

        # 3 Make person_id a one-dimensional vector
        contact_office_person_id_1D = Vector{Int64}()
        [append!(contact_office_person_id_1D, contact_office_person_id[i]) for i in 1:contact_office_N]

        # 4 Print the status of each office member
        contact_office_person_state = [df_S_[contact_office_person_id[i], :state] for i in 1:contact_office_N]

        # 5 Record the number of infected people in the office where each member works to get to the
        contact_office_I_N = [length(contact_office_person_state[i][contact_office_person_state[i] .== 2]) for i in 1:contact_office_N]
        
        # 6 Records information about each office in df_I for lambda calculation
        for i in 1:contact_office_N
            df_I_[contact_office_person_id[i], :N_pop_o] .= contact_office_I_N[i]
            df_I_[contact_office_person_id[i], :I_pop_o] .= [contact_office_person_id[i][findall(x -> x == 2, df_S_[contact_office_person_id[i], :state])] for j in 1:length(contact_office_person_id[i])]
            df_I_[contact_office_person_id[i], :N_o] .= length(contact_office_person_id[i])
        end

        list = contact_office_person_id_1D[findall(x -> x == 0, df_S_[contact_office_person_id_1D, :state])]

        contact_office_id = nothing; contact_office_person_id = nothing; contact_office_person_id_1D = nothing; 
        contact_office_person_state = nothing; contact_office_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    group_office_ = nothing; I_list_ = nothing; I_list_office = nothing; 

    return list
end

# Print out who in the company has been in contact with the infected person 
#> However, a office with 13+ has a intervention policy that only allows 1/3 of them to come to work.
function Office_Check_Constraint(I_list_, df_, df_S_, df_I_, group_office_, rngs_)

    I_list_office = I_list_[findall(x -> x != 0, df_[I_list_, :office])]
        
    if length(I_list_office) != 0
        
        # Output the office_num of offices with infected people. 
        # If the size of the company is 13 or more, only 1/3 of the employees go to work, (13*1/3 = 4 At least 4 people go to work in a large office)
        contact_office_id = unique(df_[I_list_office, :office])
        contact_office_N = length(contact_office_id)
        contact_office_size0 = [nrow(group_office_[contact_office_id[i]+1]) for i in 1:contact_office_N]
        contact_office_size = [(contact_office_size0[i] >= 13 ? Int64(round(contact_office_size0[i]*1/3)) : contact_office_size0[i]) for i in 1:contact_office_N]

        # 2 Print the person_id of the employees who work in the office where the infected person is.
        # A office with partial attendance will randomly select an employee to come to work.
        contact_office_person_id = [group_office_[contact_office_id[i]+1][sample(rngs_, 1:contact_office_size0[i], contact_office_size[i], replace=false), :person_id] for i in 1:contact_office_N]

        contact_office_person_id_1D = Vector{Int64}()
        [append!(contact_office_person_id_1D, contact_office_person_id[i]) for i in 1:contact_office_N]
        contact_office_person_state = [df_S_[contact_office_person_id[i], :state] for i in 1:contact_office_N]
        contact_office_I_N = [length(contact_office_person_state[i][contact_office_person_state[i] .== 2]) for i in 1:contact_office_N]  
        for i in 1:contact_office_N
            df_I_[contact_office_person_id[i], :N_pop_o] .= contact_office_I_N[i]
            df_I_[contact_office_person_id[i], :I_pop_o] .= [contact_office_person_id[i][findall(x -> x == 2, df_S_[contact_office_person_id[i], :state])] for j in 1:length(contact_office_person_id[i])]
            df_I_[contact_office_person_id[i], :N_o] .= length(contact_office_person_id[i])
        end
        list = contact_office_person_id_1D[findall(x -> x == 0, df_S_[contact_office_person_id_1D, :state])]

        contact_office_id = nothing; contact_office_person_id = nothing; contact_office_person_id_1D = nothing; 
        contact_office_person_state = nothing; contact_office_I_N = nothing;
        contact_office_size0 = nothing; contact_office_size = nothing;
    else
        list = Vector{Int64}()
    end
    
    group_office_ = nothing; I_list_ = nothing; I_list_office = nothing; 

    return list
end

# Print who has been in contact with an infected person in the school
function School_Check(I_list_, stage_, df_, df_S_, df_I_, group_)
    # stage_ = (1) kindergarten, (2) elementary, (3) junior, (4) high
    cols = stage_ + 3 # 4~7
        
    # 0 Print out the students among the infected
    I_list_stu = I_list_[findall(x -> x != 0, df_[I_list_, cols])]
    
    if length(I_list_stu) != 0
        # 1 Print the class_num of classes with infected people in them.
        contact_class_id = unique(df_[I_list_stu, cols])
        contact_class_N = length(contact_class_id)
        
        # 2 Print the person_id of the members of the class with the infected person in it.
        contact_class_person_id = [group_[contact_class_id[i]+1][:, :person_id] for i in 1:contact_class_N]

        # 3 Make person_id a one-dimensional vector
        contact_class_person_id_1D = Vector{Int64}()
        [append!(contact_class_person_id_1D, contact_class_person_id[i]) for i in 1:contact_class_N]

        # 4 Print out the status of each class member
        contact_class_person_state = [df_S_[contact_class_person_id[i], :state] for i in 1:contact_class_N]

        # 5 Record the number of infected people in the class attended by each member of the class.
        contact_class_I_N = [length(contact_class_person_state[i][contact_class_person_state[i] .== 2]) for i in 1:contact_class_N]

        # 6 Record information about each class in df_I for lambda calculations
        for i in 1:contact_class_N
            df_I_[contact_class_person_id[i], :N_pop_s] .= contact_class_I_N[i]
            df_I_[contact_class_person_id[i], :I_pop_s] .= [contact_class_person_id[i][findall(x -> x == 2, df_S_[contact_class_person_id[i], :state])] for j in 1:length(contact_class_person_id[i])]
            df_I_[contact_class_person_id[i], :N_s] .= length(contact_class_person_id[i])
        end

        list = contact_class_person_id_1D[findall(x -> x == 0, df_S_[contact_class_person_id_1D, :state])]

        contact_class_id = nothing; contact_class_person_id = nothing; contact_class_person_id_1D = nothing; 
        contact_class_person_state = nothing; contact_class_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    group_ = nothing; I_list_ = nothing; I_list_stu = nothing; 

    return list
end

# Print individuals who have been in contact with the infected person in the school 
# > However, a intervention policy is applied that only 1/3 attend. (high school: 2/3)
function School_Check_Constraint(I_list_, stage_, df_, df_S_, df_I_, group_, school_group_, tt_)

    cols = stage_ + 3
    gn = Int64(3 * (stage_ - 1) + tt_ % 3 + 1) # 1, 4, 7, 9 = (1) kindergarten, (4) elementary, (7) junior, (9) high
    I_list_stu = I_list_[findall(x -> x != 0, df_[I_list_, cols])]
    
    if length(I_list_stu) != 0
        tmp_list = unique(df_[I_list_stu, cols])

        # List of classes attended today
        contact_class_id = intersect(tmp_list, school_group_[gn])
        contact_class_N = length(contact_class_id)
        contact_class_person_id = [group_[contact_class_id[i]+1][:, :person_id] for i in 1:contact_class_N]
        contact_class_person_id_1D = Vector{Int64}()
        [append!(contact_class_person_id_1D, contact_class_person_id[i]) for i in 1:contact_class_N]
        contact_class_person_state = [df_S_[contact_class_person_id[i], :state] for i in 1:contact_class_N]
        contact_class_I_N = [length(contact_class_person_state[i][contact_class_person_state[i] .== 2]) for i in 1:contact_class_N]
        for i in 1:contact_class_N
            df_I_[contact_class_person_id[i], :N_pop_s] .= contact_class_I_N[i]
            df_I_[contact_class_person_id[i], :I_pop_s] .= [contact_class_person_id[i][findall(x -> x == 2, df_S_[contact_class_person_id[i], :state])] for j in 1:length(contact_class_person_id[i])]
            df_I_[contact_class_person_id[i], :N_s] .= length(contact_class_person_id[i])
        end
        list1 = contact_class_person_id_1D[findall(x -> x == 0, df_S_[contact_class_person_id_1D, :state])]

        # Class not attending school: check only teacher
        # List of classes not attending school
        shutdown_class_id = setdiff(tmp_list, school_group_[gn])
        shutdown_class_N = length(shutdown_class_id)
        shutdown_class_person_id0 = [group_[shutdown_class_id[i]+1][:, :person_id] for i in 1:shutdown_class_N]
        shutdown_class_person_age = [group_[shutdown_class_id[i]+1][:, :age] for i in 1:shutdown_class_N]
        shutdown_class_person_id  = [shutdown_class_person_id0[i][findall(shutdown_class_person_age[i] .> 18)] for i in 1:shutdown_class_N]
        shutdown_class_person_id_1D = Vector{Int64}()
        [append!(shutdown_class_person_id_1D, shutdown_class_person_id[i]) for i in 1:shutdown_class_N]
        shutdown_class_person_state = [df_S_[shutdown_class_person_id[i], :state] for i in 1:shutdown_class_N]
        shutdown_class_I_N = [length(shutdown_class_person_state[i][shutdown_class_person_state[i] .== 2]) for i in 1:shutdown_class_N]
        for i in 1:shutdown_class_N
            df_I_[shutdown_class_person_id[i], :N_pop_s] .= shutdown_class_I_N[i]
            df_I_[shutdown_class_person_id[i], :I_pop_s] .= [shutdown_class_person_id[i][findall(x -> x == 2, df_S_[shutdown_class_person_id[i], :state])] for j in 1:length(shutdown_class_person_id[i])]
            df_I_[shutdown_class_person_id[i], :N_s] .= length(shutdown_class_person_id[i])
        end
        list2 = shutdown_class_person_id_1D[findall(x -> x == 0, df_S_[shutdown_class_person_id_1D, :state])]

        list = union(list1, list2)

        contact_class_id = nothing; contact_class_person_id = nothing; contact_class_person_id_1D = nothing; 
        contact_class_person_state = nothing; contact_class_I_N = nothing;
        shutdown_class_id = nothing; shutdown_class_person_id = nothing; shutdown_class_person_id_1D = nothing; 
        shutdown_class_person_state = nothing; shutdown_class_I_N = nothing;
        tmp_list = nothing; list1 = nothing; list2 = nothing;
    else
        list = Vector{Int64}()
    end
    
    group_ = nothing; I_list_ = nothing; I_list_stu = nothing; school_group_ = nothing;

    return list
end

# Print infected people and contacts in your religion
function Religion_Check(I_list_, stage_, df_, df_S_, df_I_, df_D_, group_, rngs_)
    religion_rate = [0.8, 0.1, 0.1]

    # Religion should reflect attendance...!!!
    df_D_[:, :religion_out] .= 0

    # stage_ = (1) Christianity, (2) Catholicism, (3) Buddhism    
    cols = stage_ + 11 # 12~14 religion number # check dataframe order
    
    # 0 Print out the religious affiliation of the infected
    I_list_rel = I_list_[findall(x -> x != 0, df_[I_list_, cols])]

    if length(I_list_rel) != 0
        
        # 1 Print out the eligibility_number of the religion number of the infected person.
        contact_religion_id = unique(df_[I_list_rel, cols])
        contact_religion_N = length(contact_religion_id)
        
        # 2-0 Determine who is attending the religious institution.
        # 2-1 Print out the person_id of any member of the house of worship who attends the house of worship where the infected person is present that day, noting that the infected person may not be present.
        for i in 1:contact_religion_N
	        group_num = nrow(group_[contact_religion_id[i]+1])
            df_D_[group_[contact_religion_id[i]+1][:, :person_id], :religion_out] .= [rand(rngs_) < religion_rate[stage_] ? 1 : 0 for j in 1:group_num]
        end
        contact_religion_person_id = [group_[contact_religion_id[i]+1][findall(x->x == 1, df_D_[group_[contact_religion_id[i]+1][:, :person_id], :religion_out]), :person_id] for i in 1:contact_religion_N]

        # 3 Make person_id a one-dimensional vector
        contact_religion_person_id_1D = Vector{Int64}()
        [append!(contact_religion_person_id_1D, contact_religion_person_id[i]) for i in 1:contact_religion_N]

        # 4 Print out the status of each religion member
        contact_religion_person_state = [df_S_[contact_religion_person_id[i], :state] for i in 1:contact_religion_N]

        # 5 Enter the number of infected people in each member's religion.
        contact_religion_I_N = [length(contact_religion_person_state[i][contact_religion_person_state[i] .== 2]) for i in 1:contact_religion_N]

        # 6 Record information about each religion meeting in df_I for lambda calculation
        for i in 1:contact_religion_N
            df_I_[contact_religion_person_id[i], :N_pop_r] .= contact_religion_I_N[i]
            df_I_[contact_religion_person_id[i], :I_pop_r] .= [contact_religion_person_id[i][findall(x -> x == 2, df_S_[contact_religion_person_id[i], :state])] for j in 1:length(contact_religion_person_id[i])]
            df_I_[contact_religion_person_id[i], :N_r] .= length(contact_religion_person_id[i])
        end

        list = contact_religion_person_id_1D[findall(x -> x == 0, df_S_[contact_religion_person_id_1D, :state])]

        contact_religion_id = nothing; contact_religion_person_id = nothing; contact_religion_person_id_1D = nothing; 
        contact_religion_person_state = nothing; contact_religion_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    group_ = nothing; I_list_ = nothing; I_list_rel = nothing; 

    return list
end

# List who in your friends group has been in contact with an infected person.
function Friends_Check(I_list_, df_, df_S_, df_I_, df_D_, rngs_)

    # Determine who goes out with probability
    out_rate = 0.3
    df_D_[:, :friends_out] .= 0
    df_D_[:, :friends_m_num] .= 0
    
    # Determine who is infected and going out
    df_D_[I_list_, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(I_list_)]
    I_list_frd0 = I_list_[findall((df_D_[I_list_, :friends_out] .== 1) .& (df_[I_list_, :age] .> 2))]
    n_I = length(I_list_frd0)
    
    # Output the infected person's friends list as a one-dimensional array
    all_friends = Vector{Int64}()
    [append!(all_friends, df_[I_list_frd0[i], :friends]) for i in 1:n_I]

    # Determine if the infected person's friends are going out
    df_D_[all_friends, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(all_friends)]
    
    # Print the meeting number of the infected person and their friends
    n_max = 1
    for i in 1:n_I
        
        # If I'm out and there's no meeting yet, print the meeting number
        if (df_D_[I_list_frd0[i], :friends_out] == 1) & (df_D_[I_list_frd0[i], :friends_m_num] == 0)
        
            # Print a list of my friends who are going out today and haven't decided to meet up yet
            tmp = df_[I_list_frd0[i], :friends]
            list = tmp[findall((df_D_[tmp, :friends_out] .== 1) .& (df_D_[tmp, :friends_m_num] .== 0))]
            n_list = length(list)
            
            # How many people are in the meeting today, n_out + 1 (including me)
            n_out = Friends_Meeting_Number(n_list, rngs_)

            # Select as many as n_out from the list and put the meeting number in :friends_m_num, modifying the out status to 2 for those who decided to go out.
            # Cancel the outing if there is no one to meet
            if (n_list == 0) || (n_out == 0)
                df_D_[I_list_frd0[i], :friends_out] = 0
            else
                sample_rows = sample(rngs_, 1:n_list, n_out, replace=false)
                out_list = list[sample_rows]
                df_D_[out_list, :friends_out] .= 2
                df_D_[out_list, :friends_m_num] .= n_max
                df_D_[I_list_frd0[i], :friends_out] = 2
                df_D_[I_list_frd0[i], :friends_m_num] = n_max
                n_max = n_max + 1            
            end

            tmp = nothing;
        end
    end 
    
    group_friends = groupby(df_D_, [:friends_m_num])    

    # 0 Print out the infected people who are going out and print them again because they may cancel their outings because there is no meeting.
    I_list_frd = I_list_frd0[findall((df_D_[I_list_frd0, :friends_out] .!= 0) .& (df_D_[I_list_frd0, :friends_m_num] .!= 0))]

    if length(I_list_frd) != 0
        
        # 1 Print out the friends meeting num of friends meetings the infected person goes to.
        contact_friends_id = unique(df_D_[I_list_frd, :friends_m_num])
        contact_friends_N = length(contact_friends_id)

        # 2 Print the person_id of each friend group member who joined the infected person's friend group.
        contact_friends_person_id = [group_friends[contact_friends_id[i]+1][:, :person_id] for i in 1:contact_friends_N]
 
        # 3 Make person_id a one-dimensional vector
        contact_friends_person_id_1D = Vector{Int64}()
        [append!(contact_friends_person_id_1D, contact_friends_person_id[i]) for i in 1:contact_friends_N]

        # 4 Print the status of each friends meeting member
        contact_friends_person_state = [df_S_[contact_friends_person_id[i], :state] for i in 1:contact_friends_N]

        # 5 Record the number of infected people in each member's FRIENDS MEETING by entering the number of infected people in each member's FRIENDS MEETING.
        contact_friends_I_N = [length(contact_friends_person_state[i][contact_friends_person_state[i] .== 2]) for i in 1:contact_friends_N]
        
        # 6 Record information about each friends meeting in df_I for lambda calculation
        for i in 1:contact_friends_N
            df_I_[contact_friends_person_id[i], :N_pop_f] .= contact_friends_I_N[i]
            df_I_[contact_friends_person_id[i], :I_pop_f] .= [contact_friends_person_id[i][findall(x -> x == 2, df_S_[contact_friends_person_id[i], :state])] for j in 1:length(contact_friends_person_id[i])]
            df_I_[contact_friends_person_id[i], :N_f] .= length(contact_friends_person_id[i])
        end

        list = contact_friends_person_id_1D[findall(x -> x == 0, df_S_[contact_friends_person_id_1D, :state])]

        contact_friends_id = 1; contact_friends_person_id = nothing; contact_friends_person_id_1D = nothing; 
        contact_friends_person_state = nothing; contact_friends_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    I_list_ = nothing; I_list_fr0 = nothing; all_friends = nothing; 
    group_friends = nothing; I_list_frd = nothing;

    return list
end

# Print out who in your friends group has been in contact with an infected person 
# > However, a intervention policy is in effect that limits the size of meetings by region.
function Friends_Check_Constraint(I_list_, df_, df_S_, df_I_, df_D_, df_C_, rngs_)
    out_rate = 0.3;    df_D_[:, :friends_out] .= 0;    df_D_[:, :friends_m_num] .= 0;
    df_D_[I_list_, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(I_list_)]
    I_list_frd0 = I_list_[findall((df_D_[I_list_, :friends_out] .== 1) .& (df_[I_list_, :age] .> 2))]
    n_I = length(I_list_frd0)
    all_friends = Vector{Int64}()
    [append!(all_friends, df_[I_list_frd0[i], :friends]) for i in 1:n_I]
    df_D_[all_friends, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(all_friends)]
    
    n_max = 1
    for i in 1:n_I
        if (df_D_[I_list_frd0[i], :friends_out] == 1) & (df_D_[I_list_frd0[i], :friends_m_num] == 0)
            tmp = df_[I_list_frd0[i], :friends]
            list = tmp[findall((df_D_[tmp, :friends_out] .== 1) .& (df_D_[tmp, :friends_m_num] .== 0))]
            n_list = length(list)
            
            # Due to intervention policies, meeting size varies depending on where you live.
            n_out = Friends_Meeting_Number_Constraint(df_C_, df_[I_list_frd0[i], :region_code], n_list, rngs_)

            if (n_list == 0) || (n_out == 0)
                df_D_[I_list_frd0[i], :friends_out] = 0
            else
                sample_rows = sample(rngs_, 1:n_list, n_out, replace=false)
                out_list = list[sample_rows]
                df_D_[out_list, :friends_out] .= 2
                df_D_[out_list, :friends_m_num] .= n_max
                df_D_[I_list_frd0[i], :friends_out] = 2
                df_D_[I_list_frd0[i], :friends_m_num] = n_max
                n_max = n_max + 1            
            end
            tmp = nothing;
        end
    end 
    
    group_friends = groupby(df_D_, [:friends_m_num])    
    I_list_frd = I_list_frd0[findall((df_D_[I_list_frd0, :friends_out] .!= 0) .& (df_D_[I_list_frd0, :friends_m_num] .!= 0))]  
    if length(I_list_frd) != 0
        contact_friends_id = unique(df_D_[I_list_frd, :friends_m_num])
        contact_friends_N = length(contact_friends_id)
        contact_friends_person_id = [group_friends[contact_friends_id[i]+1][:, :person_id] for i in 1:contact_friends_N]
        contact_friends_person_id_1D = Vector{Int64}()
        [append!(contact_friends_person_id_1D, contact_friends_person_id[i]) for i in 1:contact_friends_N]
        contact_friends_person_state = [df_S_[contact_friends_person_id[i], :state] for i in 1:contact_friends_N]
        contact_friends_I_N = [length(contact_friends_person_state[i][contact_friends_person_state[i] .== 2]) for i in 1:contact_friends_N]
        for i in 1:contact_friends_N
            df_I_[contact_friends_person_id[i], :N_pop_f] .= contact_friends_I_N[i]
            df_I_[contact_friends_person_id[i], :I_pop_f] .= [contact_friends_person_id[i][findall(x -> x == 2, df_S_[contact_friends_person_id[i], :state])] for j in 1:length(contact_friends_person_id[i])]
            df_I_[contact_friends_person_id[i], :N_f] .= length(contact_friends_person_id[i])
        end
        list = contact_friends_person_id_1D[findall(x -> x == 0, df_S_[contact_friends_person_id_1D, :state])]

        contact_friends_id = 1; contact_friends_person_id = nothing; contact_friends_person_id_1D = nothing; 
        contact_friends_person_state = nothing; contact_friends_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    I_list_ = nothing; I_list_fr0 = nothing; all_friends = nothing; 
    group_friends = nothing; I_list_frd = nothing;

    return list
end

function Friends_Check_Constraint_Lv3(I_list_, df_, df_S_, df_I_, df_D_, rngs_)
    out_rate = 0.3;    df_D_[:, :friends_out] .= 0;    df_D_[:, :friends_m_num] .= 0;
    df_D_[I_list_, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(I_list_)]
    I_list_frd0 = I_list_[findall((df_D_[I_list_, :friends_out] .== 1) .& (df_[I_list_, :age] .> 2))]
    n_I = length(I_list_frd0)
    all_friends = Vector{Int64}()
    [append!(all_friends, df_[I_list_frd0[i], :friends]) for i in 1:n_I]
    df_D_[all_friends, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(all_friends)]
    
    n_max = 1
    for i in 1:n_I
        if (df_D_[I_list_frd0[i], :friends_out] == 1) & (df_D_[I_list_frd0[i], :friends_m_num] == 0)
            tmp = df_[I_list_frd0[i], :friends]
            list = tmp[findall((df_D_[tmp, :friends_out] .== 1) .& (df_D_[tmp, :friends_m_num] .== 0))]
            n_list = length(list)
            
            n_out = Friends_Meeting_Number_Constraint_Lv3(n_list, rngs_)

            if (n_list == 0) || (n_out == 0)
                df_D_[I_list_frd0[i], :friends_out] = 0
            else
                sample_rows = sample(rngs_, 1:n_list, n_out, replace=false)
                out_list = list[sample_rows]
                df_D_[out_list, :friends_out] .= 2
                df_D_[out_list, :friends_m_num] .= n_max
                df_D_[I_list_frd0[i], :friends_out] = 2
                df_D_[I_list_frd0[i], :friends_m_num] = n_max
                n_max = n_max + 1            
            end
            tmp = nothing;
        end
    end 
    group_friends = groupby(df_D_, [:friends_m_num])    
    I_list_frd = I_list_frd0[findall((df_D_[I_list_frd0, :friends_out] .!= 0) .& (df_D_[I_list_frd0, :friends_m_num] .!= 0))]  
    if length(I_list_frd) != 0
        contact_friends_id = unique(df_D_[I_list_frd, :friends_m_num])
        contact_friends_N = length(contact_friends_id)
        contact_friends_person_id = [group_friends[contact_friends_id[i]+1][:, :person_id] for i in 1:contact_friends_N]
        contact_friends_person_id_1D = Vector{Int64}()
        [append!(contact_friends_person_id_1D, contact_friends_person_id[i]) for i in 1:contact_friends_N]
        contact_friends_person_state = [df_S_[contact_friends_person_id[i], :state] for i in 1:contact_friends_N]
        contact_friends_I_N = [length(contact_friends_person_state[i][contact_friends_person_state[i] .== 2]) for i in 1:contact_friends_N]
        for i in 1:contact_friends_N
            df_I_[contact_friends_person_id[i], :N_pop_f] .= contact_friends_I_N[i]
            df_I_[contact_friends_person_id[i], :I_pop_f] .= [contact_friends_person_id[i][findall(x -> x == 2, df_S_[contact_friends_person_id[i], :state])] for j in 1:length(contact_friends_person_id[i])]
            df_I_[contact_friends_person_id[i], :N_f] .= length(contact_friends_person_id[i])
        end

        list = contact_friends_person_id_1D[findall(x -> x == 0, df_S_[contact_friends_person_id_1D, :state])]

        contact_friends_id = 1; contact_friends_person_id = nothing; contact_friends_person_id_1D = nothing; 
        contact_friends_person_state = nothing; contact_friends_I_N = nothing;
    else
        list = Vector{Int64}()
    end
    
    I_list_ = nothing; I_list_fr0 = nothing; all_friends = nothing; 
    group_friends = nothing; I_list_frd = nothing;

    return list
end

# Print who came into contact with an infected person through a chance encounter
function Encounter_Check(I_list_, df_, df_S_, df_I_, df_D_, rngs_)
    # Determine who goes out with probability
    out_rate = 0.3;    out_avg = 6;
    df_D_[:, :encounter_out] .= 0
    df_D_[:, :encounter_m_num] .= 0
    group_now_region = groupby(df_D_, [:now_region])
    
    # Choose someone to have a casual encounter with today, in an area where the infected person is today
    n_max = 1
    now_region_list = df_D_[I_list_, :now_region]
    now_region_list = unique(now_region_list)
    for i in 1:length(now_region_list)
        rc = now_region_list[i]
        n = nrow(group_now_region[rc])
        
        group_now_region[rc][:, :encounter_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:n]
        id_tmp = group_now_region[rc][:, :person_id]
        
        # People in the area who will be out and about today
        id_out = id_tmp[findall(df_D_[id_tmp, :encounter_out] .== 1)]
        n_out = length(id_out)
        
        # Number of meetings required, average meeting size is 6
        n_num = ceil(Int64, n_out / out_avg)
        
        df_D_[id_out, :encounter_m_num] .= (rand(rngs_, 1:n_num, n_out) .+ n_max)
        n_max = maximum(df_D_[id_out, :encounter_m_num])
    end
    
    group_encounter = groupby(df_D_, [:encounter_m_num])
    mc = length(group_encounter)
    for i in 2:mc
        group_encounter[i][:, :encounter_m_num] .= i-1
    end
    group_encounter = groupby(df_D_, [:encounter_m_num])

    # 0 Print out the infected people who are going out and print them again because they may cancel their outings because there is no meeting.
    I_list_en = I_list_[findall((df_D_[I_list_, :encounter_out] .!= 0) .& (df_D_[I_list_, :encounter_m_num] .!= 0))]

    if length(I_list_en) != 0
        
        # 1 Print the encounter meeting num of friend meetings the infected person goes to.
        contact_encounter_id = unique(df_D_[I_list_en, :encounter_m_num])
        contact_encounter_N = length(contact_encounter_id)

        # 2 Print the person_id of each casual meeting member who joined the infected person's casual meeting.
        contact_encounter_person_id = [group_encounter[contact_encounter_id[i]+1][:, :person_id] for i in 1:contact_encounter_N]

        # 3 Make person_id a one-dimensional vector
        contact_encounter_person_id_1D = Vector{Int64}()
        [append!(contact_encounter_person_id_1D, contact_encounter_person_id[i]) for i in 1:contact_encounter_N]

        # 4 Print the status of each encounter meeting member
        contact_encounter_person_state = [df_S_[contact_encounter_person_id[i], :state] for i in 1:contact_encounter_N]

        # 5 Enter the number of infected people in the encounter meeting that each member attends to get the
        contact_encounter_I_N = [length(contact_encounter_person_state[i][contact_encounter_person_state[i] .== 2]) for i in 1:contact_encounter_N]
        
        # 6 Record information about each friends meeting in df_I for lambda calculation
        for i in 1:contact_encounter_N
            df_I_[contact_encounter_person_id[i], :N_pop_e] .= contact_encounter_I_N[i]
            df_I_[contact_encounter_person_id[i], :I_pop_e] .= [contact_encounter_person_id[i][findall(x -> x == 2, df_S_[contact_encounter_person_id[i], :state])] for j in 1:length(contact_encounter_person_id[i])]
            df_I_[contact_encounter_person_id[i], :N_e] .= length(contact_encounter_person_id[i])
        end

        list = contact_encounter_person_id_1D[findall(x -> x == 0, df_S_[contact_encounter_person_id_1D, :state])]

        contact_encounter_id = nothing; contact_encounter_person_id = nothing; contact_encounter_person_id_1D = nothing; 
        contact_encounter_person_state = nothing; contact_encounter_I_N = nothing;
    else
        list = Vector{Int64}()
    end

    I_list_ = nothing; I_list_en = nothing; group_now_region = nothing; group_encounter = nothing;
    now_region_list = nothing; id_tmp = nothing; id_out = nothing; 

    return list
end

function Compute_Lambda(id_, vs_, df_, df_S_, df_I_, beta0, rngs_)
    # Household, School/Office, Friends Gatherings, Religion, Random Encounter
    beta = [beta0, beta0, beta0, beta0, beta0, beta0]

    # Household (1) + School (2) + Workplace (3) + Religion (4) + Friends' mettion (5) + Random encounter (6) : 6 total
    lambda = zeros(Float64, 6);     n_I = zeros(Int64, 6);    vs_I = zeros(Float64, 6)
    
    n_I[1] = df_I_[id_, 3*1]
    tmp = df_I_[id_, 3*1-1]
    # The reason we do vs_count -1 is because we incremented vs_count by 1 while moving the infected to the next state first (Check_State() first).
    # Different people have different relative infectiousness.
    vs_I[1] = sum( Float64[(vs_[df_S_[tmp[i], :vs_count]-1] / df_S_[tmp[i], :vs_N]) * df_S_[tmp[i], :RI] for i in 1:length(tmp)] )
    lambda[1] = beta[1] * vs_I[1] / df_I_[id_, 3*1+1]^0.8
    for i in 2:3
        n_I[i] = df_I_[id_, 3*i]
        tmp = df_I_[id_, 3*i-1]
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_[id_, i+15]
    end
    for i in 4:6
        n_I[i] = df_I_[id_, 3*i]
        tmp = df_I_[id_, 3*i-1]
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_I_[id_, 3*i+1]
    end 
    sum_lambda = NaNMath.sum(lambda)

    if sum_lambda > rand(rngs_)
        df_S_[id_, :I_type] = findall(n_I .!= 0)
        n = 1
    else 
        n = 0
    end
    
    lambda = nothing; n_I = nothing; vs_I = nothing; beta = nothing;

    return n
end

# National capital area Private Gathering Size Limits, beta_5/2.0 (Intervention policy Lv.2)
function Compute_Lambda_Constraint_Lv2(id_, vs_, df_, df_S_, df_I_, beta0, df_C_, rngs_)
    beta = [beta0, beta0, beta0, beta0, beta0, beta0]

    lambda = zeros(Float64, 6);     n_I = zeros(Int64, 6);    vs_I = zeros(Float64, 6)
    
    # Household
    n_I[1] = df_I_[id_, 3*1];    tmp = df_I_[id_, 3*1-1];
    vs_I[1] = sum( Float64[(vs_[df_S_[tmp[i], :vs_count]-1] / df_S_[tmp[i], :vs_N]) * df_S_[tmp[i], :RI] for i in 1:length(tmp)] )
    lambda[1] = beta[1] * vs_I[1] / df_I_[id_, 3*1+1]^0.8
    # School/Workplace
    for i in 2:3
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_[id_, i+15]
    end
    # Religion/Random Encounter
    for i in [4, 6]
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_I_[id_, 3*i+1]
    end 
    # Friends' Meeting
    n_I[5] = df_I_[id_, 3*5];    tmp = df_I_[id_, 3*5-1];
    vs_I[5] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
    lambda[5] = df_C_[df_[id_, :region_code], :beta] * beta[5] * vs_I[5] / df_I_[id_, 3*5+1]

    sum_lambda = NaNMath.sum(lambda)

    if sum_lambda > rand(rngs_)
        df_S_[id_, :I_type] = findall(n_I .!= 0)
        n = 1
    else 
        n = 0
    end
    
    lambda = nothing; n_I = nothing; vs_I = nothing; beta = nothing;

    return n
end

# Entire Country Private Gathering Size Limits, beta_5/2.5 (Intervention policy Lv.3)
function Compute_Lambda_Constraint_Lv3(id_, vs_, df_, df_S_, df_I_, beta0, rngs_)
    beta = [beta0, beta0, beta0, beta0, beta0/2.5, beta0]
    lambda = zeros(Float64, 6);     n_I = zeros(Int64, 6);    vs_I = zeros(Float64, 6)

    n_I[1] = df_I_[id_, 3*1];    tmp = df_I_[id_, 3*1-1];
    vs_I[1] = sum( Float64[(vs_[df_S_[tmp[i], :vs_count]-1] / df_S_[tmp[i], :vs_N]) * df_S_[tmp[i], :RI] for i in 1:length(tmp)] )
    lambda[1] = beta[1] * vs_I[1] / df_I_[id_, 3*1+1]^0.8
    for i in 2:3
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_[id_, i+15]
    end
    for i in 4:6
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_I_[id_, 3*i+1]
    end 
    sum_lambda = NaNMath.sum(lambda)

    if sum_lambda > rand(rngs_)
        df_S_[id_, :I_type] = findall(n_I .!= 0)
        n = 1
    else 
        n = 0
    end
    
    lambda = nothing; n_I = nothing; vs_I = nothing; beta = nothing;

    return n
end

function Check_State(id_, t_, df_S_)
    s = df_S_[id_, :state]
    if df_S_[id_, s+7] == df_S_[id_, :day_count]
        if s == 2
            df_S_[id_, :date_I] = t_ - df_S_[id_, :period_I]
            df_S_[id_, :date_R] = t_
            df_S_[id_, :vs_count] = df_S_[id_, :vs_count] + 1
        end
        df_S_[id_, :day_count] = 1
        df_S_[id_, :state] = df_S_[id_, :state] + 1
    else
        df_S_[id_, :day_count] = df_S_[id_, :day_count] + 1
        if s == 2 # As long as you can infect others, your viral shedding will change daily.
            df_S_[id_, :vs_count] = df_S_[id_, :vs_count] + 1
        end
    end

end

# Setting up a new infected person
function Setting_I_State(I_list_, t_, df_S_, df_I_, std_, vs_A_, rngs_)
    df_S_[I_list_, :date_E] .= t_
    df_S_[I_list_, :day_count] .= 1
    df_S_[I_list_, :state] .= 1
    df_S_[I_list_, :period_E] .= [State_E_Period(rngs_) for i in 1:length(I_list_)]
    df_S_[I_list_, :period_I] .= [State_I_Period() for i in 1:length(I_list_)]
    df_S_[I_list_, :RI] .= [Relative_Infectiousness(std_, rngs_) for i in 1:length(I_list_)]

    for j in 1:length(I_list_)
        id = I_list_[j]
        
        # Not infect others (in E)
        tmp = df_S_[id, :period_E] - 2
        df_S_[id, :vs_count] = (tmp >= 1 ? 1 : (tmp == 0 ? 2 : 3))
        df_S_[id, :date_unt] = (tmp >= 1 ? tmp : 1)
        
        # Duration of infecting others
        df_S_[id, :date_inf] = 10 - df_S_[id, :vs_count] + 1
        
        # viral shedding normalization
        df_S_[id, :vs_N] = (df_S_[id, :date_inf] == 10 ? vs_A_[1] : (df_S_[id, :date_inf] == 9 ? vs_A_[2] : vs_A_[3]))

        # Record who is infected
        contact_I_list = [df_I_[id, 2+3*(i-1)] for i in 1:6]
        contact_I_list_1d = Vector{Int64}()
        [append!(contact_I_list_1d, contact_I_list[i]) for i in 1:6]
        df_S_[id, :I_pop] = contact_I_list_1d
    end

    I_list_ = nothing; contact_I_list = nothing; contact_I_list_1d = nothing;
end

# After a day, reset today's contact data frame
function Reset_I_Data(list_, df_I_)
    for i in 1:6
        j = 2+3*(i-1)
        df_I_[list_, j] .= [Vector{Int64}()]
        df_I_[list_, j+1] .= 0
        df_I_[list_, j+2] .= 0
    end
end

function InitialSetting(df_, df_S_, std_, vs_A_, rngs_)  

    dfr = DataFrame(CSV.File("./src/Initial_Setting_R_201031.csv"))  
    dfi = DataFrame(CSV.File("./src/Initial_Setting_inf_201027.csv"))
    
    A = groupby(df_, [:region_code, :age_group])

    R_list = Vector{Int64}();    I_list = Vector{Int64}();
    for rn in 1:250
        for ag in 1:18
            r, c = size(df_[(df_.region_code .== rn) .& (df_.age_group .== ag), :])
            if r != 0
                tmp = A[(region_code = rn, age_group = ag,)]    
                sample_rows = sample(rngs_, 1:nrow(tmp), dfr[rn, ag]+dfi[rn, ag], replace=false)
                append!(R_list, tmp[sample_rows[1:dfr[rn, ag]], :person_id])
                append!(I_list, tmp[sample_rows[dfr[rn, ag]+1:dfr[rn, ag]+dfi[rn, ag]], :person_id])
                tmp = 1
            end # constraint
        end # age group
    end # region

    R0 = length(R_list)
    df_S_[R_list, :state] .= 3; 
    df_S_[R_list, :period_E] .= -1;    df_S_[R_list, :period_I] .= -1;
    df_S_[R_list, :date_E] .= -1;    df_S_[R_list, :date_I] .= -1;    df_S_[R_list, :date_R] .= -1;
    df_S_[R_list, :vs_count] .= -1;    df_S_[R_list, :date_unt] .= -1;    df_S_[R_list, :date_inf] .= -1;    df_S_[R_list, :day_count] .= 1;
    df_S_[R_list, :I_type] .= [zeros(Int64, 1) for i in 1:R0];    df_S_[R_list, :I_pop] .= [zeros(Int64, 1) for i in 1:R0];
    df_S_[R_list, :RI] .= -1;
    df_S_[R_list, :vs_N] .= -1;

    I0 = length(I_list)
    df_S_[I_list, :state] .= 1;
    df_S_[I_list, :date_E] .= -1;    df_S_[I_list, :day_count] .= 1;
    df_S_[I_list, :period_E] .= [State_E_Period(rngs_) for i in 1:I0];    df_S_[I_list, :period_I] .= [State_I_Period() for i in 1:I0];
    df_S_[I_list, :day_count] .= 1;
    df_S_[I_list, :I_type] .= [zeros(Int64, 1) for i in 1:I0];    df_S_[I_list, :I_pop] .= [zeros(Int64, 1) for i in 1:I0];
    df_S_[I_list, :RI] .= [Relative_Infectiousness(std_, rngs_) for i in 1:I0];
    for i in 1:I0
        id = I_list[i]
        tmp = df_S_[id, :period_E] - 2
        df_S_[id, :vs_count] = (tmp >= 1 ? 1 : (tmp == 0 ? 2 : 3))
        df_S_[id, :date_unt] = (tmp >= 1 ? tmp : 1)
        df_S_[id, :date_inf] = 10 - df_S_[id, :vs_count] + 1
        df_S_[id, :vs_N] = (df_S_[id, :date_inf] == 10 ? vs_A_[1] : (df_S_[id, :date_inf] == 9 ? vs_A_[2] : vs_A_[3]))
    end

    dfr = nothing; dfi = nothing;
    R_list = nothing; I_list = nothing; A = nothing; sample_rows = nothing;

    GC.gc()

end

function InitialZeroSetting(df_D_, df_S_, df_I_)
    df_D_[:, 2:7] .= 0
    df_S_[:, 12] .= [Vector{Int64}()]
    df_S_[:, 13] .= [Vector{Int64}()]
    df_S_[:, 2:11] .= 0
    df_S_[:, 14:15] .= 0
    for i in 1:6
        j = 2+3*(i-1)
        df_I_[:, j] .= [Vector{Int64}()]
        df_I_[:, j+1] .= 0
        df_I_[:, j+2] .= 0
    end
end

function main()
    Random.seed!(1)
    thread_size = Threads.nthreads()
    println(" intervention policy + relative infectiousness + thread size = ", thread_size)

    df = ReadCSVFile()
    df = InfodfSetting(df)
    N, c = size(df)

    group_house = groupby(df, [:house_id]);    group_region = groupby(df, [:region_code]);
    group_pre = groupby(df, [:pre]);    group_elementary = groupby(df, [:elementary]);    group_junior = groupby(df, [:junior]);    group_high = groupby(df, [:high]);
    group_office = groupby(df, [:office]);
    group_christian = groupby(df, [:christian]);    group_catholic = groupby(df, [:catholic]);    group_buddhism = groupby(df, [:buddhism]);

    vs = InitialSetting_vs()
    vs_A = zeros(3); vs_A[1] = sum(vs[1:10]); vs_A[2] = sum(vs[2:10]); vs_A[3] = sum(vs[3:10]);
    df_C = MakeRegionConstraint()

    df_D = [MakedfDay(N) for i in 1:thread_size]
    df_S = [MakedfState(N) for i in 1:thread_size]
    df_I = [MakedfLambda(N) for i in 1:thread_size]

    ncycles = 25;    ntries1 = 38;    ntries2 = 54;    ntries3 = 65;    ntries4 = 100;    dt = 1;
    rngs = [MersenneTwister(i) for i in 1:thread_size];   
    t = ones(Int64, thread_size);    ts = ones(Int64, thread_size);    nt = Int64(thread_size*ncycles);
    I_list = [Vector{Int64}() for i in 1:thread_size];    N_I = zeros(Int64, thread_size);    ii = zeros(Int64, thread_size)
    contact_list = [Vector{Int64}() for i in 1:thread_size];    check_list = [Vector{Int64}() for i in 1:thread_size];
    other_region_list = [[] for i in 1:thread_size];    school_group = [[] for i in 1:thread_size];
    start0 = zeros(thread_size);    end0 = zeros(thread_size);
    b_rank = [0.8, 0.8, 0.8, 0.8, 0.8];
    s_rank = [0.5, 0.5, 0.5, 0.5, 0.5];

    println(" multi threads start! : total cycle = ", ncycles, " final time = ", ntries4, " beta = ", b_rank, " RI std = ", s_rank)
    Threads.@threads for rank = 1:thread_size

        for n in 1:thread_size:nt
            fn = n + rank - 1

            rs = Int64(round(time()+fn*(fn+1)))
            Random.seed!(rngs[rank], rs)

            t[rank] = 1;    ii[rank] = 0;    ts[rank] = 1;    other_region_list[rank] = [];    school_group[rank] = [];
            InitialZeroSetting(df_D[rank], df_S[rank], df_I[rank])
            other_region_list[rank] = MovePopulation(group_region, rngs[rank])
            school_group[rank] = Setting_Group_Education_Rotation(df, group_pre, group_elementary, group_junior, group_high, rngs[rank])
            
            InitialSetting(df, df_S[rank], s_rank[rank], vs_A, rngs[rank])
            start0[rank] = time()
                        
            for t[rank] in 1:ntries1-1        
                I_list[rank] = df_S[rank][findall((df_S[rank].state .== 2) .& (df_S[rank].RI != 0)), :person_id]
                N_I[rank] = length(I_list[rank])
                
                # Compute the lambda by checking only those who have been in contact with the infected person.
                contact_list[rank] = Vector{Int64}()
                
                # Household, Everyday
                append!(contact_list[rank], House_Check(I_list[rank], df, df_S[rank], df_I[rank], group_house))

                # School/Workplace, Weekday
                if 2 <= t[rank] % 7 <= 6
                    append!(contact_list[rank], Office_Check(I_list[rank], df, df_S[rank], df_I[rank], group_office))
                    append!(contact_list[rank], School_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], group_pre))
                    append!(contact_list[rank], School_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], group_elementary))
                    append!(contact_list[rank], School_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], group_junior))
                    append!(contact_list[rank], School_Check(I_list[rank], 4, df, df_S[rank], df_I[rank], group_high))
                end

                # Religion, Sunday
                if t[rank] % 7 == 1
                    append!(contact_list[rank], Religion_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], df_D[rank], group_christian, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], df_D[rank], group_catholic, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], df_D[rank], group_buddhism, rngs[rank]))
                end

                # Friends' meeting, Everyday
                append!(contact_list[rank], Friends_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))
                
                # Random Encounter, Everyday
                for ii[rank] in 1:250
                    df_D[rank][group_region[ii[rank]][:, :person_id], :now_region] .= Random.shuffle!(rngs[rank], other_region_list[rank][ii[rank]])
                end
                append!(contact_list[rank], Encounter_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))

                # Transition E and I to the next state
                check_list[rank] = df_S[rank][findall((df_S[rank].state .== 1) .|| (df_S[rank].state .== 2)), :person_id]
                [Check_State(check_list[rank][ii[rank]], t[rank], df_S[rank]) for ii[rank] in 1:length(check_list[rank])]
                
                # Find new infected person, need to remove duplicates from contact list (unique)
                contact_list[rank] = unique(contact_list[rank])
                check_list[rank] = [Compute_Lambda(contact_list[rank][ii[rank]], vs, df, df_S[rank], df_I[rank], b_rank[rank], rngs[rank]) for ii[rank] in 1:length(contact_list[rank])]
                Setting_I_State(contact_list[rank][findall(check_list[rank] .== 1)], t[rank], df_S[rank], df_I[rank], s_rank[rank], vs_A, rngs[rank])

                # Reset today's contact data for the next day
                Reset_I_Data(I_list[rank], df_I[rank]);
                Reset_I_Data(contact_list[rank], df_I[rank])
            end
            
            println(" rank = ", rank, " ", fn, " : Intervention Policy Level 1 --------------------- ")

            for t[rank] in ntries1:ntries2-1

                I_list[rank] = df_S[rank][findall((df_S[rank].state .== 2) .& (df_S[rank].RI != 0)), :person_id]
                N_I[rank] = length(I_list[rank])
                contact_list[rank] = Vector{Int64}()
                
                append!(contact_list[rank], House_Check(I_list[rank], df, df_S[rank], df_I[rank], group_house))

                if 2 <= t[rank] % 7 <= 6
                    append!(contact_list[rank], Office_Check_Constraint(I_list[rank], df, df_S[rank], df_I[rank], group_office, rngs[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 1, df, df_S[rank], df_I[rank], group_pre, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 2, df, df_S[rank], df_I[rank], group_elementary, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 3, df, df_S[rank], df_I[rank], group_junior, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 4, df, df_S[rank], df_I[rank], group_high, school_group[rank], ts[rank]))
                    ts[rank] += 1
                end
                
                if t[rank] % 7 == 1
                    append!(contact_list[rank], Religion_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], df_D[rank], group_christian, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], df_D[rank], group_catholic, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], df_D[rank], group_buddhism, rngs[rank]))
                end

                append!(contact_list[rank], Friends_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))
                
                for ii[rank] in 1:250
                    df_D[rank][group_region[ii[rank]][:, :person_id], :now_region] .= Random.shuffle!(rngs[rank], other_region_list[rank][ii[rank]])
                end
                append!(contact_list[rank], Encounter_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))

                check_list[rank] = df_S[rank][findall((df_S[rank].state .== 1) .|| (df_S[rank].state .== 2)), :person_id]
                [Check_State(check_list[rank][ii[rank]], t[rank], df_S[rank]) for ii[rank] in 1:length(check_list[rank])]
                
                contact_list[rank] = unique(contact_list[rank])
                check_list[rank] = [Compute_Lambda(contact_list[rank][ii[rank]], vs, df, df_S[rank], df_I[rank], b_rank[rank], rngs[rank]) for ii[rank] in 1:length(contact_list[rank])]
                Setting_I_State(contact_list[rank][findall(check_list[rank] .== 1)], t[rank], df_S[rank], df_I[rank], s_rank[rank], vs_A, rngs[rank])
            
                Reset_I_Data(I_list[rank], df_I[rank])
                Reset_I_Data(contact_list[rank], df_I[rank])
                
            end
            
            println(" rank = ", rank, " ", fn, " : Intervention Policy Level 2 --------------------- ")

            for t[rank] in ntries2:ntries3-1

                I_list[rank] = df_S[rank][findall((df_S[rank].state .== 2) .& (df_S[rank].RI != 0)), :person_id]
                N_I[rank] = length(I_list[rank])
                
                contact_list[rank] = Vector{Int64}()
                
                append!(contact_list[rank], House_Check(I_list[rank], df, df_S[rank], df_I[rank], group_house))

                if 2 <= t[rank] % 7 <= 6
                    append!(contact_list[rank], Office_Check_Constraint(I_list[rank], df, df_S[rank], df_I[rank], group_office, rngs[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 1, df, df_S[rank], df_I[rank], group_pre, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 2, df, df_S[rank], df_I[rank], group_elementary, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 3, df, df_S[rank], df_I[rank], group_junior, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 4, df, df_S[rank], df_I[rank], group_high, school_group[rank], ts[rank]))
                    ts[rank] += 1
                end
                
                if t[rank] % 7 == 1
                    append!(contact_list[rank], Religion_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], df_D[rank], group_christian, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], df_D[rank], group_catholic, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], df_D[rank], group_buddhism, rngs[rank]))
                end

                append!(contact_list[rank], Friends_Check_Constraint(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], df_C, rngs[rank]))
                
                for ii[rank] in 1:250
                    df_D[rank][group_region[ii[rank]][:, :person_id], :now_region] .= Random.shuffle!(rngs[rank], other_region_list[rank][ii[rank]])
                end
                append!(contact_list[rank], Encounter_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))

                check_list[rank] = df_S[rank][findall((df_S[rank].state .== 1) .|| (df_S[rank].state .== 2)), :person_id]
                [Check_State(check_list[rank][ii[rank]], t[rank], df_S[rank]) for ii[rank] in 1:length(check_list[rank])]
                
                contact_list[rank] = unique(contact_list[rank])
                check_list[rank] = [Compute_Lambda_Constraint_Lv2(contact_list[rank][ii[rank]], vs, df, df_S[rank], df_I[rank], b_rank[rank], df_C, rngs[rank]) for ii[rank] in 1:length(contact_list[rank])]
                Setting_I_State(contact_list[rank][findall(check_list[rank] .== 1)], t[rank], df_S[rank], df_I[rank], s_rank[rank], vs_A, rngs[rank])
            
                Reset_I_Data(I_list[rank], df_I[rank])
                Reset_I_Data(contact_list[rank], df_I[rank])                
            end

            println(" rank = ", rank, " ", fn, " : Intervention Policy Level 3 --------------------- ")

            for t[rank] in ntries3:ntries4

                I_list[rank] = df_S[rank][findall((df_S[rank].state .== 2) .& (df_S[rank].RI != 0)), :person_id]
                N_I[rank] = length(I_list[rank])
                
                contact_list[rank] = Vector{Int64}()
                
                append!(contact_list[rank], House_Check(I_list[rank], df, df_S[rank], df_I[rank], group_house))

                if 2 <= t[rank] % 7 <= 6
                    append!(contact_list[rank], Office_Check_Constraint(I_list[rank], df, df_S[rank], df_I[rank], group_office, rngs[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 1, df, df_S[rank], df_I[rank], group_pre, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 2, df, df_S[rank], df_I[rank], group_elementary, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 3, df, df_S[rank], df_I[rank], group_junior, school_group[rank], ts[rank]))
                    append!(contact_list[rank], School_Check_Constraint(I_list[rank], 4, df, df_S[rank], df_I[rank], group_high, school_group[rank], ts[rank]))
                    ts[rank] += 1
                end
                
                if t[rank] % 7 == 1
                    append!(contact_list[rank], Religion_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], df_D[rank], group_christian, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], df_D[rank], group_catholic, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], df_D[rank], group_buddhism, rngs[rank]))
                end

                append!(contact_list[rank], Friends_Check_Constraint_Lv3(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))
                
                for ii[rank] in 1:250
                    df_D[rank][group_region[ii[rank]][:, :person_id], :now_region] .= Random.shuffle!(rngs[rank], other_region_list[rank][ii[rank]])
                end
                append!(contact_list[rank], Encounter_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))

                check_list[rank] = df_S[rank][findall((df_S[rank].state .== 1) .|| (df_S[rank].state .== 2)), :person_id]
                [Check_State(check_list[rank][ii[rank]], t[rank], df_S[rank]) for ii[rank] in 1:length(check_list[rank])]
                
                contact_list[rank] = unique(contact_list[rank])
                check_list[rank] = [Compute_Lambda_Constraint_Lv3(contact_list[rank][ii[rank]], vs, df, df_S[rank], df_I[rank], b_rank[rank], rngs[rank]) for ii[rank] in 1:length(contact_list[rank])]
                Setting_I_State(contact_list[rank][findall(check_list[rank] .== 1)], t[rank], df_S[rank], df_I[rank], s_rank[rank], vs_A, rngs[rank])
            
                Reset_I_Data(I_list[rank], df_I[rank])
                Reset_I_Data(contact_list[rank], df_I[rank])
                
            end


            end0[rank] = time()

            println(" rank = ", rank, " ", fn, " : OneSet complete! >>>>>>>>> total time taken = ", end0[rank]-start0[rank], " s")
            CSV.write("./IBM_ivp_b$(b_rank[rank])_std$(s_rank[rank])_$(rs).csv", df_S[rank][(df_S[rank].state .!= 0), :])

            GC.gc()
        end
        
    end
    
    t = nothing; N_I = nothing; ii = nothing; b_rank = nothing; s_rank = nothing; ts = nothing;
    start0 = nothing; end0 = nothing;
    I_list = nothing; contact_list = nothing; check_list = nothing; 
    other_region_list = nothing; school_group = nothing;
    df = nothing; df_S = nothing; df_I = nothing; df_D = nothing; df_C = nothing;
    group_house = nothing; group_region = nothing;
    group_pre = nothing; group_elementary = nothing; group_junior = nothing; group_high = nothing;
    group_office = nothing; group_christian = nothing; group_catholic = nothing; group_buddhism = nothing;
    GC.gc()

end

main()
