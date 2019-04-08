using LinearAlgebra
using LinQuadOptInterface
const LQOI = LinQuadOptInterface
using MathOptInterface
const MOI = MathOptInterface

const status_dict = Dict(
    :Unbounded => MOI.DUAL_INFEASIBLE,
    :Infeasible => MOI.INFEASIBLE,
    :Optimal => MOI.OPTIMAL,
    :Unknown => MOI.OTHER_ERROR
)

const p_status_dict = Dict(
    :Unbounded => MOI.NO_SOLUTION,
    :Infeasible => MOI.NO_SOLUTION,
    :Optimal => MOI.FEASIBLE_POINT,
    :Unknown => MOI.NO_SOLUTION
)

const con_dict = Dict(
      Cchar('G') => '>',
      Cchar('L') => '<',
      Cchar('E') => '=',
)

function Optimizer(instance)

    #=
        Read relevant data
    =#
    params = LQOI.get_parameters(instance)
    sense = LQOI.get_optimization_sense(instance)
    c = LQOI.get_objective_coefficients(instance)
    Qobj = LQOI.get_objective_quadratic_coefficients(instance)
    @assert isempty(Qobj)
    Qcon = LQOI.get_constraint_quadratic_coefficients(instance)
    for Q in Qcon
        @assert isempty(Q)
    end
    A = LQOI.get_constraint_coefficients(instance)
    b = LQOI.get_constraint_constant(instance)
    # LQOI.get_constraint_constant_range_ub(instance)
    contype = LQOI.get_constraint_type(instance)
    lb = LQOI.get_variable_lower_bound(instance)
    ub = LQOI.get_variable_upper_bound(instance)
    vartype = LQOI.get_variable_type(instance)
    for v in vartype
        @assert v == Cchar('C')
    end
    sos = LQOI.get_sos_constraints(instance)
    @assert length(sos) == 0

    #=
        Modify data for the internal format
    =#
    new_c = vcat(rationalize.(c), -rationalize.(c))
    rat_A = rationalize.(A)
    new_A = hcat(rat_A, -rat_A)

    real_ub = zeros(Int, 0)
    real_lb = zeros(Int, 0)
    for i in 1:length(ub)
        if ub[i] < Inf
                push!(real_ub, i)
        end
        if lb[i] > -Inf
                push!(real_lb, i)
        end
    end
    Id = Matrix{Rational{Int}}(I, length(ub), length(ub))

    new_A = hvcat((1, 2, 2), 
        new_A, 
        Id[real_ub,:], -Id[real_ub,:], 
        Id[real_lb,:], -Id[real_lb,:])

    new_b = vcat(rationalize.(b), 
                rationalize.(ub[real_ub]),
                rationalize.(lb[real_lb]))

    new_contype = vcat(
        map(x -> con_dict[x], contype),
        fill('<', length(real_ub)),
        fill('>', length(real_lb))
    )

    # for i in 1:length(new_b)
    #     if new_b[i] < 0.0
    #         new_b[i] .*= -1
    #         if new_contype[i] == '<'
    #             new_contype[i] = '>'
    #         elseif new_contype[i] == '>'
    #             new_contype[i] = '<'
    #         end
    #         new_A[i,:] .*= -1
    #     end
    # end

    new_sense = sense == :minimize ? :Min : :Max

    #=
        Optimize!
    =#
    status, x = simplex(new_c, new_sense, new_A, new_b, new_contype)

    #=
    Send reults back to JuMP and MOI
    =#
    xx = x[1:div(end,2)] + x[1+div(end,2):end]
    solution = LQOI.DataLinQuadOptimizerSolution(
        status_dict[status],
        p_status_dict[status],
        MOI.NO_SOLUTION,
        xx,
        fill(NaN, length(xx)),
        A*xx,
        fill(NaN, length(b)),
        zeros(0),
        zeros(0),
        dot(c, xx),
        NaN,
        1,
        1,
        0
    )

    return solution
end