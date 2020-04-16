module HopfieldRBM
struct RBM{T}
    N::Integer
    P::Integer
    ξ::AbstractMatrix{T}
    m::AbstractVector{T}
    σ::AbstractVector{T}
    function RBM{T}(N,P) where T
        ξ=ones(T,N,P)
        m=ones(T,P)
        σ=ones(T,N)
        new{T}(N,P,ξ,m,σ)
    end
    RBM(N,P)=RBM{Float64}(N,P)
    function RBM(ξ::AbstractMatrix{T}) where T
        N,P=size(ξ)
        m=similar(ξ,P)
        σ=similar(ξ,N)
        σ.=1
        S=new{T}(N,P,ξ,m,σ)
        refresh_m!(S)
        S
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

function glauber_flip_σ!(S::RBM{T}) where T
    H=energy(S)
    for i in 1:S.N
        flip_σ!(S,i)
        fH=energy(S)
        p=1 / (1 + exp(H - fH))
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
        glauber_flip_σ!(S)
        for μ in 1:S.P, k in 1:S.N
            g[k,μ]-=S.m[μ]*S.σ[k]/NC
        end
    end
end
using PyPlot

function testfigure()
    S=RBM(rand((-1.0,1.0),100,2).*0.01)
    pastξ=copy(S.ξ)
    s1=ones(100);
    s1[1:50].*=-1;
    s2=copy(s1);
    s2[1:25].*=-1;
    s2[75:100].*=-1;
    g=zeros(100,2);

    for  i in 1:40
        g.=0
        ∇logP!(S,g,s1)
        ∇logP!(S,g,s2)
        S.ξ .+= 0.05 .* g
    end
    S.ξ

    pygui(true)
    subplot(1,2,1)
    imshow(pastξ, interpolation="nearest", aspect="auto")
    subplot(1,2,2)
    imshow(S.ξ, interpolation="nearest", aspect="auto")
end

testfigure()
end # module
