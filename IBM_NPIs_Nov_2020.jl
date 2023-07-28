""" 
    2023-08-01
    IBM Julia language code developed by Min-Kyung Chae. (in NIMS)
    multi-thread code
    a.out > julia --threads 2 filename.jl

    individual-based model (IBM) to simulate the spread of COVID-19 in South Korea

    intervention policy : November 2020 - February 2022
    Level 1 : school capacity reduce + large workplace telework
    Level 2 : private gathering size is maximum 4 in the national capital area
    Level 3 : private gathering size is maximum 4 in the entire country
"""

using Base.Threads
using DataFrames, CSV, StatsBase, Statistics, Distributions
using Random
using Plots, Images
using RData
using NaNMath
import CodecBzip2


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
    out_max = 9  # 모임 규는 +1 을 해야함. 나 포함. (즉, 최대 규모 10)
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Friends_Meeting_Number_Constraint(df_C_, rc_, n_, rngs_)
    out_max = df_C_[rc_, :constraint] - 1  # 모임 규는 +1 을 해야함. 나 포함.
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Friends_Meeting_Number_Constraint_Lv3(n_, rngs_)
    out_max = 3  # 모임 규는 +1 을 해야함. 나 포함. (즉, 최대 규모 4)
    max = (n_ < out_max ? n_ : out_max)
    rn = (max == 0 ? 0 : rand(rngs_, 1:max))
    
    return rn
end

function Setting_Group_Education_Rotation(df_, group1_, group2_, group3_, group4_, rngs_)
    # 각 지역의 학급 중에 30%만 등교를 해 (1/3 씩 등교) >  등교 그룹을 지정해 주는 거임
    # 250 개 시군구를 들려야 해서 초콤 오래 걸림, 1번만 하면 되는 거라 그냥 두고 있음

    k = 1;     school_group_ = [];
    for cols in 4:7 # 유, 초, 중, 고 
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

# 모임 제한 데이터 프레임
function MakeRegionConstraint()
    df_C_ = DataFrame("region_code" => 1:250, "constraint" => 10, "beta" => 1.0) 
    # 수도권만 4인
    df_C_[1:25, :constraint] .= 4; df_C_[1:25, :beta] .= 0.5;  # 서울
    df_C_[50:59, :constraint] .= 4; df_C_[50:59, :beta] .= 0.5;  # 인천
    df_C_[75:116, :constraint] .= 4; df_C_[75:116, :beta] .= 0.5; # 경기

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

# 가구 내 감염자와 접촉한 사람 출력
function House_Check(I_list_, df_, df_S_, df_I_, group_house_)

    # 1 감염된 사람이 있는 가구의 house_id 를 출력해
    contact_house_id = unique(df_[I_list_, :house_id])
    contact_house_N = length(contact_house_id)
    
    # 2 감염된 사람이 있는 가구에 살고 있는 각 가구 구성원들의 person_id 를 출력해
    contact_house_person_id = [group_house_[contact_house_id[i]][:, :person_id] for i in 1:contact_house_N]
    
    # 3 person_id 를 1차원 vector로 만들어 줌 (나중에 contact_list 추가를 편하게 하기 위해서)
    contact_house_person_id_1D = Vector{Int64}()
    [append!(contact_house_person_id_1D, contact_house_person_id[i]) for i in 1:contact_house_N]
    
    # 4 각 가구 구성원의 상태를 출력해
    contact_house_person_state = [df_S_[contact_house_person_id[i], :state] for i in 1:contact_house_N]
    
    # 5 각 구성원들이 살고 있는 가구에서 감염된 사람의 수를 입력해
    contact_house_I_N = [length(contact_house_person_state[i][contact_house_person_state[i] .== 2]) for i in 1:contact_house_N]
    
    # 6 lambda 계산을 위해 df_I 에 각 person 관한 정보를 기록함
    # N_pop_h : 가정내 감염자 수 , I_pop_h : 가정내 감염자의 person id , N_h : 가정내 인구 수
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

# 회사 내 감염자와 접촉한 사람 출력
function Office_Check(I_list_, df_, df_S_, df_I_, group_office_)

    # 0 감염된 사람중 회사원을 출력해
    I_list_office = I_list_[findall(x -> x != 0, df_[I_list_, :office])]
    
    if length(I_list_office) != 0
        # 1 감염된 사람이 있는 회사의 office_num 를 출력해
        contact_office_id = unique(df_[I_list_office, :office])
        contact_office_N = length(contact_office_id)

        # 2 감염된 사람이 있는 office를 다니는 있는 구성원들의 person_id 를 출력해
        contact_office_person_id = [group_office_[contact_office_id[i]+1][:, :person_id] for i in 1:contact_office_N]

        # 3 person_id 를 1차원 vector로 만들어줌
        contact_office_person_id_1D = Vector{Int64}()
        [append!(contact_office_person_id_1D, contact_office_person_id[i]) for i in 1:contact_office_N]

        # 4 각 office 구성원의 상태를 출력해
        contact_office_person_state = [df_S_[contact_office_person_id[i], :state] for i in 1:contact_office_N]

        # 5 각 구성원들이 다니는 office에서 감염된 사람의 수를 입력해
        contact_office_I_N = [length(contact_office_person_state[i][contact_office_person_state[i] .== 2]) for i in 1:contact_office_N]
        
        # 6 lambda 계산을 위해 df_I에 각 office 관한 정보를 기록함
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

# 회사 내 감염자와 접촉한 사람 출력 > 그런데, 13명 이사의 회사의 경우 1/3 만 출근하는 방역 정책 적용.
function Office_Check_Constraint(I_list_, df_, df_S_, df_I_, group_office_, rngs_)

    I_list_office = I_list_[findall(x -> x != 0, df_[I_list_, :office])]
        
    if length(I_list_office) != 0
        
        # 1 감염된 사람이 있는 회사의 office_num 를 출력. 
        # 회사 규모가 13이상이면 1/3 만 출근 하고, (13*1/3 = 4 큰회사는 최소 4명이 출근쓰..)
        contact_office_id = unique(df_[I_list_office, :office])
        contact_office_N = length(contact_office_id)
        contact_office_size0 = [nrow(group_office_[contact_office_id[i]+1]) for i in 1:contact_office_N]
        contact_office_size = [(contact_office_size0[i] >= 13 ? Int64(round(contact_office_size0[i]*1/3)) : contact_office_size0[i]) for i in 1:contact_office_N]

        # 2 감염된 사람이 있는 office를 다니는 있는 구성원들의 person_id 를 출력해
        # 일부만 출근 하는 회사는 출근할 사원을 랜덤으로 선택해 준다
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

# 학교 내 감염자와 접촉한 사람 출력
function School_Check(I_list_, stage_, df_, df_S_, df_I_, group_)
    # stage_ = (1) kindergarten, (2) elementary, (3) junior, (4) high
    cols = stage_ + 3 # 4~7
        
    # 0 감염된 사람중 학생을 출력해
    I_list_stu = I_list_[findall(x -> x != 0, df_[I_list_, cols])]
    
    if length(I_list_stu) != 0
        # 1 감염된 사람이 있는 class의 class_num 를 출력해
        contact_class_id = unique(df_[I_list_stu, cols])
        contact_class_N = length(contact_class_id)
        
        # 2 감염된 사람이 있는 class를 다니는 있는 구성원들의 person_id 를 출력해
        contact_class_person_id = [group_[contact_class_id[i]+1][:, :person_id] for i in 1:contact_class_N]

        # 3 person_id 를 1차원 vector로 만들어줌
        contact_class_person_id_1D = Vector{Int64}()
        [append!(contact_class_person_id_1D, contact_class_person_id[i]) for i in 1:contact_class_N]

        # 4 각 class 구성원의 상태를 출력해
        contact_class_person_state = [df_S_[contact_class_person_id[i], :state] for i in 1:contact_class_N]

        # 5 각 구성원들이 다니는 class에서 감염된 사람의 수를 입력해
        contact_class_I_N = [length(contact_class_person_state[i][contact_class_person_state[i] .== 2]) for i in 1:contact_class_N]

        # 6 lambda 계산을 위해 df_I에 각 class 관한 정보를 기록함
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

# 학교 내 감염자와 접촉한 사람 출력 > 그런데, 유초중 1/3, 고 2/3 만 등교하는 방역 정책 적용.
function School_Check_Constraint(I_list_, stage_, df_, df_S_, df_I_, group_, school_group_, tt_)

    cols = stage_ + 3
    gn = Int64(3 * (stage_ - 1) + tt_ % 3 + 1) # 1, 4, 7, 9 = (1) kindergarten, (4) elementary, (7) junior, (9) high
    I_list_stu = I_list_[findall(x -> x != 0, df_[I_list_, cols])]
    
    if length(I_list_stu) != 0
        # 1 감염된 사람이 있는 class의 class_num 를 출력해
        tmp_list = unique(df_[I_list_stu, cols])

        # 오늘 등교한 반
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

        # 등교 안한 반 : only teacher
        # 등교를 하지 않는 반의 list
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

# 종교내 감염자와 첩촉한 사람 출력
function Religion_Check(I_list_, stage_, df_, df_S_, df_I_, df_D_, group_, rngs_)
    religion_rate = [0.8, 0.1, 0.1]

    # 얘네는 종교에 참석율을 반영해 줘야함..!!!
    df_D_[:, :religion_out] .= 0

    # stage_ = (1) Christianity, (2) Catholicism, (3) Buddhism    
    cols = stage_ + 11 # 12~14 religion number (df 를 합치면서 순서가 바뀜, 확인해야함.)
    
    # 0 감염된 사람중 종교인을 출력해
    I_list_rel = I_list_[findall(x -> x != 0, df_[I_list_, cols])]

    if length(I_list_rel) != 0
        
        # 1 감염된 사람이 있는 종교 number의 religion_number 를 출력해
        contact_religion_id = unique(df_[I_list_rel, cols])
        contact_religion_N = length(contact_religion_id)
        
        # 2-0 해당 종교 시설에서 참석할 사람을 결정해
        # 2-1 감염된 사람이 있는 종교 시설을 다니는 있는 구성원들 중에 그날 출석한 person_id 를 출력해, 감염자가 참석 안 할 수도 있음.
        for i in 1:contact_religion_N
	        group_num = nrow(group_[contact_religion_id[i]+1])
            df_D_[group_[contact_religion_id[i]+1][:, :person_id], :religion_out] .= [rand(rngs_) < religion_rate[stage_] ? 1 : 0 for j in 1:group_num]
        end
        contact_religion_person_id = [group_[contact_religion_id[i]+1][findall(x->x == 1, df_D_[group_[contact_religion_id[i]+1][:, :person_id], :religion_out]), :person_id] for i in 1:contact_religion_N]

        # 3 person_id 를 1차원 vector로 만들어줌
        contact_religion_person_id_1D = Vector{Int64}()
        [append!(contact_religion_person_id_1D, contact_religion_person_id[i]) for i in 1:contact_religion_N]

        # 4 각 religion 구성원의 상태를 출력해
        contact_religion_person_state = [df_S_[contact_religion_person_id[i], :state] for i in 1:contact_religion_N]

        # 5 각 구성원들이 다니는 religion 감염된 사람의 수를 입력해
        contact_religion_I_N = [length(contact_religion_person_state[i][contact_religion_person_state[i] .== 2]) for i in 1:contact_religion_N]

        # 6 lambda 계산을 위해 df_I에 각 religion meeting 관한 정보를 기록함
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

# 친구 모임내 감염자와 접촉한 사람 출력
function Friends_Check(I_list_, df_, df_S_, df_I_, df_D_, rngs_)
    out_rate = 0.3

    # 확률로 외출 의사 반영해 줘야함..!!!
    df_D_[:, :friends_out] .= 0
    df_D_[:, :friends_m_num] .= 0
    
    # 감염자 중에 외출하는 사람 결정
    df_D_[I_list_, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(I_list_)]
    I_list_frd0 = I_list_[findall((df_D_[I_list_, :friends_out] .== 1) .& (df_[I_list_, :age] .> 2))]
    n_I = length(I_list_frd0)
    
    # 감염자의 친구 목록 1차원 배열로 출력
    all_friends = Vector{Int64}()
    [append!(all_friends, df_[I_list_frd0[i], :friends]) for i in 1:n_I]

    # 감염자의 친구가 외출하는지 결정
    df_D_[all_friends, :friends_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:length(all_friends)]
    
    # 감염자와 감염자 친구들의 모임 번호를 출력
    n_max = 1
    for i in 1:n_I
        
        # 내가 외출하는데, 아직 모임이 없으면 모임 번호 출력
        if (df_D_[I_list_frd0[i], :friends_out] == 1) & (df_D_[I_list_frd0[i], :friends_m_num] == 0)
        
            # 오늘 외출하고, 아직 모임이 결정되지 않은 내 친구 목록 출력
            tmp = df_[I_list_frd0[i], :friends]
            list = tmp[findall((df_D_[tmp, :friends_out] .== 1) .& (df_D_[tmp, :friends_m_num] .== 0))]
            n_list = length(list)
            
            # 오늘 모임은 몇 명인가, n_out + 1 (나 포함)
            n_out = Friends_Meeting_Number(n_list, rngs_)

            # list 에서 n_out 만큼 선택 해서 :friends_m_num에 모임 번호 넣기, 외출이 결정된 애들은 외출 여부를 2로 수정
            # 만날 사람이 없으면 외출 취소
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

    # 0 감염된 사람중 외출 하는 사람을 출력해, 다시 출력하는 이유는 모임이 없어서 외출을 취소 할 수 도 있어서임
    I_list_frd = I_list_frd0[findall((df_D_[I_list_frd0, :friends_out] .!= 0) .& (df_D_[I_list_frd0, :friends_m_num] .!= 0))]

    if length(I_list_frd) != 0
        
        # 1 감염된 사람이 가는 친구 모임의 friends meeting num 를 출력해
        contact_friends_id = unique(df_D_[I_list_frd, :friends_m_num])
        contact_friends_N = length(contact_friends_id)

        # 2 감염된 사람이 가는 친구 모임에 참가한 각 친구 모임 구성원들의 person_id 를 출력해
        contact_friends_person_id = [group_friends[contact_friends_id[i]+1][:, :person_id] for i in 1:contact_friends_N]
 
        # 3 person_id 를 1차원 vector로 만들어줌
        contact_friends_person_id_1D = Vector{Int64}()
        [append!(contact_friends_person_id_1D, contact_friends_person_id[i]) for i in 1:contact_friends_N]

        # 4 각 friends meeting 구성원의 상태를 출력해
        contact_friends_person_state = [df_S_[contact_friends_person_id[i], :state] for i in 1:contact_friends_N]

        # 5 각 구성원들이 다니는 friends meeting에서 감염된 사람의 수를 입력해
        contact_friends_I_N = [length(contact_friends_person_state[i][contact_friends_person_state[i] .== 2]) for i in 1:contact_friends_N]
        
        # 6 lambda 계산을 위해 df_I에 각 friends meeting 관한 정보를 기록함
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

# 친구 모임내 감염자와 접촉한 사람 출력 > 그런데, 지역 별로 만남 규모에 제한이 있는 방역 정책 효과 적용.
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
            
            # 방역 정책 때문에, 거주 지역마다 모임 규모가 달라짐.
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
            
            # 오늘 모임은 몇 명인가, n_out + 1 (나 포함)
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

# 우연한 만남으로 감염자와 접촉한 사람 출력
function Encounter_Check(I_list_, df_, df_S_, df_I_, df_D_, rngs_)
    out_rate = 0.3;    out_avg = 6;

    df_D_[:, :encounter_out] .= 0
    df_D_[:, :encounter_m_num] .= 0
    group_now_region = groupby(df_D_, [:now_region])
    
    # 감염자가 오늘 있는 지역에서, 오늘 우연한 만남을 할 사람 선택  
    n_max = 1
    now_region_list = df_D_[I_list_, :now_region]
    now_region_list = unique(now_region_list)
    for i in 1:length(now_region_list)
        rc = now_region_list[i]
        n = nrow(group_now_region[rc])
        
        group_now_region[rc][:, :encounter_out] .= [rand(rngs_) < out_rate ? 1 : 0 for i in 1:n]
        id_tmp = group_now_region[rc][:, :person_id]
        
        # 그 지역에 오늘 외출할 사람들
        id_out = id_tmp[findall(df_D_[id_tmp, :encounter_out] .== 1)]
        n_out = length(id_out)
        
        # 필요한 모임의 수, 평균 모임 크기는 6
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

    # 0 감염된 사람중 외출 하는 사람을 출력해, 다시 출력하는 이유는 모임이 없어서 외출을 취소 할 수 도 있어서임
    I_list_en = I_list_[findall((df_D_[I_list_, :encounter_out] .!= 0) .& (df_D_[I_list_, :encounter_m_num] .!= 0))]

    if length(I_list_en) != 0
        
        # 1 감염된 사람이 가는 친구 모임의 encounter meeting num 를 출력해
        contact_encounter_id = unique(df_D_[I_list_en, :encounter_m_num])
        contact_encounter_N = length(contact_encounter_id)

        # 2 감염된 사람이 가는 우연한 모임에 참가한 각 우연한 모임 구성원들의 person_id 를 출력해
        contact_encounter_person_id = [group_encounter[contact_encounter_id[i]+1][:, :person_id] for i in 1:contact_encounter_N]

        # 3 person_id 를 1차원 vector로 만들어줌
        contact_encounter_person_id_1D = Vector{Int64}()
        [append!(contact_encounter_person_id_1D, contact_encounter_person_id[i]) for i in 1:contact_encounter_N]

        # 4 각 encounter meeting 구성원의 상태를 출력해
        contact_encounter_person_state = [df_S_[contact_encounter_person_id[i], :state] for i in 1:contact_encounter_N]

        # 5 각 구성원들이 다니는 encounter meeting에서 감염된 사람의 수를 입력해
        contact_encounter_I_N = [length(contact_encounter_person_state[i][contact_encounter_person_state[i] .== 2]) for i in 1:contact_encounter_N]
        
        # 6 lambda 계산을 위해 df_I에 각 friends meeting 관한 정보를 기록함
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
    beta = [beta0, beta0, beta0, beta0, beta0, beta0] # 가구, 학교/회사, 친구모임, 종교, 우연한모임

    # 가구(1) + 학교(2) + 직장(3) + 종교(4) + 친구(5) + 우연(6) : 총 6개
    lambda = zeros(Float64, 6);     n_I = zeros(Int64, 6);    vs_I = zeros(Float64, 6)
    
    n_I[1] = df_I_[id_, 3*1]
    tmp = df_I_[id_, 3*1-1]
    # 여기서 vs_count -1 을 하는 이유는 감염자들을 먼저 다음 state 로 옮기면서 vs_count 를 1씩 증가 시켰기 때문.  (Check_State 를 먼저 함)
    # 사람마다 relative infectiousness 가 다름
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

# 수도권 사적 모임 인원 제한, beta/2.0
function Compute_Lambda_Constraint_Lv2(id_, vs_, df_, df_S_, df_I_, beta0, df_C_, rngs_)
    beta = [beta0, beta0, beta0, beta0, beta0, beta0]

    # 가구(1) + 학교(2) + 직장(3) + 종교(4) + 친구(5) + 우연(6) : 총 6개
    lambda = zeros(Float64, 6);     n_I = zeros(Int64, 6);    vs_I = zeros(Float64, 6)
    
    # 집
    n_I[1] = df_I_[id_, 3*1];    tmp = df_I_[id_, 3*1-1];
    vs_I[1] = sum( Float64[(vs_[df_S_[tmp[i], :vs_count]-1] / df_S_[tmp[i], :vs_N]) * df_S_[tmp[i], :RI] for i in 1:length(tmp)] )
    lambda[1] = beta[1] * vs_I[1] / df_I_[id_, 3*1+1]^0.8
    # 학교/회사
    for i in 2:3
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_[id_, i+15]
    end
    #종교, 우연한 만남
    for i in [4, 6]
        n_I[i] = df_I_[id_, 3*i];        tmp = df_I_[id_, 3*i-1];
        vs_I[i] = sum( Float64[(vs_[df_S_[tmp[j], :vs_count]-1] / df_S_[tmp[j], :vs_N]) * df_S_[tmp[j], :RI] for j in 1:length(tmp)] )
        lambda[i] = beta[i] * vs_I[i] / df_I_[id_, 3*i+1]
    end 
    #친구
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

# 전국 사적 모임 인원 제한, beta_5/2.5 
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
        if s == 2 # 남을 감염 시킬 수 있는 상태면, viral shedding 도 매일 바뀜.
            df_S_[id_, :vs_count] = df_S_[id_, :vs_count] + 1
        end
    end

end

# 새로운 감염자 셋팅
function Setting_I_State(I_list_, t_, df_S_, df_I_, std_, vs_A_, rngs_)
    df_S_[I_list_, :date_E] .= t_
    df_S_[I_list_, :day_count] .= 1
    df_S_[I_list_, :state] .= 1
    df_S_[I_list_, :period_E] .= [State_E_Period(rngs_) for i in 1:length(I_list_)]
    df_S_[I_list_, :period_I] .= [State_I_Period() for i in 1:length(I_list_)]
    df_S_[I_list_, :RI] .= [Relative_Infectiousness(std_, rngs_) for i in 1:length(I_list_)]

    for j in 1:length(I_list_)
        id = I_list_[j]
        
        # E 기간에서 남을 감염시키지 않는 기간
        tmp = df_S_[id, :period_E] - 2
        df_S_[id, :vs_count] = (tmp >= 1 ? 1 : (tmp == 0 ? 2 : 3))
        df_S_[id, :date_unt] = (tmp >= 1 ? tmp : 1)
        
        # 남을 감염 시키는 기간
        df_S_[id, :date_inf] = 10 - df_S_[id, :vs_count] + 1
        
        # viral shedding normalization
        df_S_[id, :vs_N] = (df_S_[id, :date_inf] == 10 ? vs_A_[1] : (df_S_[id, :date_inf] == 9 ? vs_A_[2] : vs_A_[3]))

        # 누구에게 감염 되었는지를 기록
        contact_I_list = [df_I_[id, 2+3*(i-1)] for i in 1:6]
        contact_I_list_1d = Vector{Int64}()
        [append!(contact_I_list_1d, contact_I_list[i]) for i in 1:6]
        df_S_[id, :I_pop] = contact_I_list_1d
    end

    I_list_ = nothing; contact_I_list = nothing; contact_I_list_1d = nothing;
end

# 하루 지나면, 오늘 접촉자 데이터 리셋
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
            end #조건
        end #나이 그룹
    end #지역

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
                
                # 감염된 사람과 접촉한 사람만 따로 체크해서 lambda 를 계산한다.
                contact_list[rank] = Vector{Int64}()
                
                # 가구, 매일
                append!(contact_list[rank], House_Check(I_list[rank], df, df_S[rank], df_I[rank], group_house))

                # 학교 + 회사, 월화수목금
                if 2 <= t[rank] % 7 <= 6
                    append!(contact_list[rank], Office_Check(I_list[rank], df, df_S[rank], df_I[rank], group_office))
                    append!(contact_list[rank], School_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], group_pre))
                    append!(contact_list[rank], School_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], group_elementary))
                    append!(contact_list[rank], School_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], group_junior))
                    append!(contact_list[rank], School_Check(I_list[rank], 4, df, df_S[rank], df_I[rank], group_high))
                end

                # 종교, 일요일
                if t[rank] % 7 == 1
                    append!(contact_list[rank], Religion_Check(I_list[rank], 1, df, df_S[rank], df_I[rank], df_D[rank], group_christian, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 2, df, df_S[rank], df_I[rank], df_D[rank], group_catholic, rngs[rank]))
                    append!(contact_list[rank], Religion_Check(I_list[rank], 3, df, df_S[rank], df_I[rank], df_D[rank], group_buddhism, rngs[rank]))
                end

                # 친구, 매일
                append!(contact_list[rank], Friends_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))
                
                # 우연한 만남, 매일
                for ii[rank] in 1:250
                    df_D[rank][group_region[ii[rank]][:, :person_id], :now_region] .= Random.shuffle!(rngs[rank], other_region_list[rank][ii[rank]])
                end
                append!(contact_list[rank], Encounter_Check(I_list[rank], df, df_S[rank], df_I[rank], df_D[rank], rngs[rank]))

                # E 와 I 를 다음 state로 전이
                check_list[rank] = df_S[rank][findall((df_S[rank].state .== 1) .|| (df_S[rank].state .== 2)), :person_id]
                [Check_State(check_list[rank][ii[rank]], t[rank], df_S[rank]) for ii[rank] in 1:length(check_list[rank])]
                
                # 새로운 감염자 찾기, contact list 에서 중복 제거 해야함(unique)
                contact_list[rank] = unique(contact_list[rank])
                check_list[rank] = [Compute_Lambda(contact_list[rank][ii[rank]], vs, df, df_S[rank], df_I[rank], b_rank[rank], rngs[rank]) for ii[rank] in 1:length(contact_list[rank])]
                Setting_I_State(contact_list[rank][findall(check_list[rank] .== 1)], t[rank], df_S[rank], df_I[rank], s_rank[rank], vs_A, rngs[rank])

                # 다음 날을 위해 오늘자 접촉 데이터 리셋
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
            CSV.write("./IBM3_ivp_b$(b_rank[rank])_std$(s_rank[rank])_$(rs).csv", df_S[rank][(df_S[rank].state .!= 0), :])

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
