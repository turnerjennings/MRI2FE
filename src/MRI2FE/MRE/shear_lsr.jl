module OptModel
using JuMP, Ipopt

function optimize_model(gp::Float64, gpp::Float64, w::Float64)
    """
    Calculate the prony series constants for an equivalent complex shear modulus using least squares optimization

    Arguments:
        gp (float): Storage modulus
        gpp (float): Loss modulus
        w (float): Modulus frequency

    Returns
        Ginf (float): Long-time shear modulus
        G1 (float): short-time shear modulus
        tau (float): time constant
    """

    m = Model(Ipopt.Optimizer)
    set_optimizer_attribute(m, MOI.Silent(), true)

    @variable(m, g[1:2] ≥ 0)
    @variable(m, τ ≥ 0)

    gpex = @expression(m, g[1] + ((g[2] - g[1]) * w^2 * τ^2) / (1 + w^2 * τ^2))
    gppex = @expression(m, ((g[2] - g[1]) * w * τ) / (1 + w^2 * τ^2))

    @constraint(m, gppex / gpex ≥ gpp / gp)
    @constraint(m, gppex / gpex ≤ gpp / gp)

    #@constraint(m, τ == 1/w)
    @objective(m, Min, (gp - gpex)^2 + (gpp - gppex)^2)

    optimize!(m)



    return value(g[1]), value(g[2]), value(τ)
end

end