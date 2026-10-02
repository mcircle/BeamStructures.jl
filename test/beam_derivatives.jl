@testset "Beam RHS derivatives overwrite reusable buffers" begin
    p = BS.SciMLBase.NullParameters()
    for T in (Float32, Float64)
        state = T[0.2, 0.4, 0, 0, 0.1, -0.2, 0.3]
        reference = ForwardDiff.jacobian(x -> BS.ode(x, p, zero(T)), state)
        buffer = fill(T(NaN), 7, 7)
        BS.jac!(buffer, state, p, zero(T))
        @test buffer ≈ reference
        @test BS.jac(state, p, zero(T)) ≈ reference
        lambda = T[1, 2, 3, 4, 5, 6, 7]
        vjp = fill(T(NaN), 7)
        BS.vjp_beam!(vjp, lambda, state, zero(T))
        @test vjp ≈ -reference' * lambda
        augmented = hcat(vcat(lambda, state), vcat(2lambda, state))
        matrix_vjp = fill(T(NaN), 14, 2)
        BS.vjp!(matrix_vjp, augmented, p, zero(T))
        @test matrix_vjp[1:7,1] ≈ vjp
        @test matrix_vjp[1:7,2] ≈ 2vjp
        @test matrix_vjp[8:14,1] ≈ BS.ode(state, p, zero(T))
        @test matrix_vjp[8:14,2] ≈ matrix_vjp[8:14,1]
        curved = fill(T(NaN), 6)
        BS.vjp_curved_beam!(curved, lambda[1:6], state, zero(T))
        @test curved ≈ -reference[1:6,1:6]' * lambda[1:6]
    end
end

@testset "Empty force orientation preserves beam tangent" begin
    beam = BS.Beam(25.0,1.0,5.0,0.0)
    beams = (;Beam_1=beam)
    tangent = CRC.zero_tangent(beams)
    state = zeros(7,2,1)
    state_tangent = zero(state)
    force_tangent = zeros(3)
    @test BS.forcesbackatend!(state_tangent,tangent,force_tangent,state,beams,Int[]) === tangent
    @test BS.forcesbackatstart!(state_tangent,tangent,force_tangent,state,beams,Int[]) === tangent
end

function one_sided_branch_model(parameters,orientation)
    beam(offset) = BS.Beam(parameters[offset:offset+6]...)
    beams = (Beam_1=beam(1),Beam_2=beam(8),Beam_3=beam(15))
    if orientation === :incoming
        nodes = (Node_1=BS.Clamp(0.0,0.0,0.0,0.0,0.0,0.0),
                 Node_2=BS.Clamp(50.0,0.0,0.0,0.0,0.0,0.0),
                 Node_3=BS.Branch(25.0,25.0,0.0,0.0,0.0,0.0),
                 Node_4=BS.Clamp(50.0,50.0,0.0,0.0,0.0,0.0))
        structure = BS.Structure([0.0 0 1 1;0 0 1 0;1 1 0 0;1 0 0 0])
    else
        nodes = (Node_1=BS.Clamp(0.0,0.0,0.0,0.0,0.0,0.0),
                 Node_2=BS.Branch(25.0,25.0,0.0,0.0,0.0,0.0),
                 Node_3=BS.Clamp(50.0,0.0,0.0,0.0,0.0,0.0),
                 Node_4=BS.Clamp(50.0,50.0,0.0,0.0,0.0,0.0))
        structure = BS.Structure([0.0 0 1 0;0 0 1 1;1 1 0 0;0 1 0 0])
    end
    structure,beams,nodes
end

function one_sided_branch_loss(parameters,orientation)
    structure,beams,nodes = one_sided_branch_model(parameters,orientation)
    states = reshape(parameters[22:end],7,2,3)
    residual = zeros(eltype(parameters),3,4)
    weights = reshape(collect(1.0:12.0),3,4)
    sum(BS.residuals!(residual,structure,states,beams,nodes) .* weights)
end

function one_sided_pipeline_loss(parameters,orientation)
    structure,beams,nodes = one_sided_branch_model(parameters,orientation)
    initial = reshape(parameters[22:end],3,4)
    solutions,_,_ = structure(initial,beams,nodes)
    sum(solutions)
end

@testset "One-sided Branch pullback agrees with ForwardDiff" begin
    beams = (BS.Beam(35.0,1.0,5.0,0.01;θs=0.2),
             BS.Beam(45.0,1.0,5.0,-0.01;θs=-0.1),
             BS.Beam(30.0,1.0,5.0,0.005;θs=0.3))
    parameters = vcat(collect.(Tuple.(beams))...,collect(range(-0.2,0.3;length=42)))
    for orientation in (:incoming,:outgoing)
        objective = parameters -> one_sided_branch_loss(parameters,orientation)
        reverse = only(Zygote.gradient(objective,parameters))
        forward = ForwardDiff.gradient(objective,parameters)
        @info "one-sided Branch gradient comparison" orientation beam_max_abs=maximum(abs,reverse[1:21]-forward[1:21]) state_max_abs=maximum(abs,reverse[22:end]-forward[22:end])
        @test reverse ≈ forward rtol=1e-8 atol=1e-10
    end
end


@testset "Full Structure pipeline agrees with ForwardDiff" begin
    beams = (BS.Beam(35.0,1.0,5.0,0.01;θs=0.2),
             BS.Beam(45.0,1.0,5.0,-0.01;θs=-0.1),
             BS.Beam(30.0,1.0,5.0,0.005;θs=0.3))
    initial = 0.01 .* randn(MersenneTwister(20261002),12)
    parameters = vcat(collect.(Tuple.(beams))...,initial)
    for orientation in (:incoming,:outgoing)
        objective = parameters -> one_sided_pipeline_loss(parameters,orientation)
        reverse = only(Zygote.gradient(objective,parameters))
        forward = ForwardDiff.gradient(objective,parameters)
        beam_difference = abs.(reverse[1:21] .- forward[1:21])
        beam_index = argmax(beam_difference)
        @info "full Structure pipeline gradient comparison" orientation beam_max_abs=beam_difference[beam_index] beam=cld(beam_index,7) field=fieldnames(typeof(beams[1]))[mod1(beam_index,7)] reverse_value=reverse[beam_index] forward_value=forward[beam_index] state_max_abs=maximum(abs,reverse[22:end]-forward[22:end])
        @test all(isapprox.(reverse,forward;rtol=2e-3,atol=1e-5))
    end
end
