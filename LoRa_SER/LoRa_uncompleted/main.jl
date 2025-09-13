using Distributed
@everywhere using Plots
@everywhere using FFTW
@everywhere using ProgressMeter
@everywhere using Statistics
@everywhere using LinearAlgebra


#julia -p auto main.jl 
@everywhere function nom_dist(mean::Float64, var::Float64)#正規分布を返す乱数(複素数)，::Float64は書いても書かなくてもパフォーマンスには影響しない
    x=rand()
    y=rand()
    Re=sqrt(-2.0*var*log(x))*cos(2pi*y)+mean
    Im=sqrt(-2.0*var*log(x))*sin(2pi*y)+mean
    return Re+im*Im
end
@everywhere function add_noise(mean::Float64,var::Float64)
    Noise=nom_dist(mean,var/2.0)
    return Noise
end
@everywhere function dft(y::Array,M::Int64;ovs=1)
    Y=zeros(ComplexF64,M)#周波数領域
    for n in 0:M-1
        for k in 0:M-1 #配列は1から始まることに注意
            Y[n+1] += y[k+1]*exp(-im*2pi*n*k/M)
        end
    end
    return Y
end

@everywhere function demodu(sf::Int64,bw::Int64,tc::Float64,ts::Float64,x::Array,Y::Array,Y_abs::Array,d::Array)
    M=2^sf
    for k in 1:M
        d[k]=x[k]*exp(-im*2pi*((k-1)*tc)^2*bw/(2*ts))
    end
    Y=dft(d,M)
    Y_abs=1/(M)*abs.(Y)
    m_star=findfirst(x->x==maximum(Y_abs[1:2^sf]),Y_abs[1:2^sf])-1
    return m_star
end

@everywhere function ser_calc(snr::Float64,sf::Int64,bw::Int64,iter::Int64)
    M=2^sf #number of chips
    tc=1/bw # chip length [sec]
    ts=tc*M # symbol length [sec]
    var = 1.0/10^(snr/10.0)
    error::Int64=0
    x=zeros(ComplexF64,M) # fscc (LoRa) signal on equivalent low-pass system
    Y=zeros(ComplexF64,M)
    Y_abs=zeros(ComplexF64,M)
    d=zeros(ComplexF64,M)
    m=0
    m_star=0
    for i in 1:iter
       m=rand(0:2^sf-1) 
        for k in 1:M
            x[k]=exp(im*2pi*mod((m+(k-1))*tc,M*tc)^2*bw/(2*ts))
        end
        for k in 1:M
            x[k]+=add_noise(0.0,var)
        end
        m_star=demodu(sf,bw,tc,ts,x,Y,Y_abs,d)
        if m!=m_star
            error+=1
        end
    end
    GC.gc()
    return error
end



function main(sf::Int64,bw::Int64,resolution::Float64,iter::Int64)
    if !isdir("result")#ディレクトリ生成
        mkdir("result")
    end
    cd("./result")#ディレクトリ移動
    minSNR=-20.0
    maxSNR=0.0
    amount_plot=Int64((maxSNR-minSNR)/resolution+1)
    multi=500
    println("The parameters of the program are as follows... Spreading Factor=$sf  , Band Width=$bw")
    ser_array= zeros(Float64,amount_plot,2)##snr,serが入る
    j = 1#カウンタ
    snr=minSNR
    while snr <= maxSNR
        _iter=iter
        _iter=Int64(_iter/multi)
        A=@showprogress pmap(x->ser_calc(snr,sf,bw,_iter),1:multi)
        println("$snr done")
        error = 0
        for x in A
            error += x
        end
        ser_array[j,1] =snr
        ser_array[j,2] = Float64(error)/iter
        j+=1
        snr+=resolution
    end
    io1 = open("single_ser_sf$(sf)_minSNR$(minSNR)_maxSNR$(maxSNR)_resolution$(resolution)_iter$(iter).csv", "w")
    for i in 1 : amount_plot
        println(io1, ser_array[i,1], "," , ser_array[i,2])
    end
    close(io1)
    #axis_x=range(minSNR,maxSNR,length=amount_plot)
    
    # p1=plot(ser_array[:,1],ser_array[:,2],
    #     legend=false,
    #     xgrid=false,
    #     ygrid=false,
    #     fmt=:png)
    # savefig(p1,"result_ser_sf$(sf)_minSNR$(minSNR)_maxSNR$(maxSNR)_resolution$(resolution)_iter$(iter).png")
    # display(p1)
    # println("Press Enter to exit this program (or function)")
    # readline(stdin)
end
@elapsed main(9,125*10^3,0.5,10^3)