module HopfieldRBM
struct RBM{T}
    N::Integer
    P::Integer
    ξ::AbstractMatrix{T}
    m::AbstractVector{T}
    σ::AbstractVector{T}
    function RBM{T}(N::Integer,P::Integer) where T
        ξ=ones(T,N,P)
        m=ones(T,P)
        σ=ones(T,N)
        new{T}(N,P,ξ,m,σ)
    end
    RBM(N::Integer,P::Integer)=RBM{Float64}(N,P)
    function RBM(ξ::AbstractMatrix{T},σ::AbstractVector{T}) where T
        N,P=size(ξ)
        N==length(σ) || error("σ must be long $N, not ",length(σ),".")
        m=similar(σ,P)
        S=new{T}(N,P,ξ,m,σ)
        refresh_m!(S)
        S
    end
    function RBM(ξ::AbstractMatrix{T}) where T
        N,P=size(ξ)
        σ=similar(ξ,N)
        σ.=1
        RBM(ξ,σ)
    end
end

function refresh_m!(S::RBM{T}) where T
    S.m .= 0
    for μ in 1:S.P
        for i in 1:S.N
            S.m[μ]+=sign(S.ξ[i,μ]*S.σ[i])
        end
        S.m[μ]/=S.N
    end
end

function energy(S::RBM{T}) where T
    S.N * sum(abs2, S.m) / 2
end

function apply_σ!(S::RBM{T},σ) where T
    S.σ .= σ
    refresh_m!(S)
end

function flip_σ!(S::RBM{T},k) where T
    for μ in 1:S.P
        S.m[μ]=round(S.m[μ]*S.N-2sign(S.ξ[k,μ]*S.σ[k]))/S.N
    end
    S.σ[k]=-sign(S.σ[k])
end

function glauber_flip_σ!(S::RBM{T},β::T) where T
    H=energy(S)
    for i in 1:S.N
        flip_σ!(S,i)
        fH=energy(S)
        p=1 / (1 + exp(β*(H - fH)))
        if rand() < p
            H=fH
        else
            flip_σ!(S,i)
        end
    end
end

function ∇logP!(S::RBM{T},g::AbstractMatrix{T},σD::AbstractVector{T}) where T
    apply_σ!(S,σD)
    for μ in 1:S.P, k in 1:S.N
        g[k,μ]+=S.m[μ]*S.σ[k]
    end
    NC=5
    for i in 1:NC
        glauber_flip_σ!(S,one(T))
        for μ in 1:S.P, k in 1:S.N
            g[k,μ]-=S.m[μ]*S.σ[k]/NC
        end
    end
end
using PyPlot

function testfigure()
    S=RBM(rand((-1.0,1.0),100,3).*0.1)
    pastξ=copy(S.ξ)
    s1=ones(100);
    s1[1:50].*=-1;
    s2=copy(s1);
    s2[1:25].*=-1;
    s2[75:100].*=-1;
    s3=copy(s2);
    s3[1:12].*=-1;
    s3[37:50].*=-1;
    s3[51:100].=s3[1:50]
    g=similar(S.ξ);

    for  i in 1:50
        g.=0
        ∇logP!(S,g,s1)
        ∇logP!(S,g,s2)
        ∇logP!(S,g,s3)
        S.ξ .+= 0.05 .* g
    end

    pygui(true)
    subplot(1,2,1)
    imshow(pastξ, interpolation="nearest", aspect="auto")
    colorbar()
    subplot(1,2,2)
    imshow(S.ξ, interpolation="nearest", aspect="auto")
    colorbar()
    tight_layout()
end

testfigure()
end # module
