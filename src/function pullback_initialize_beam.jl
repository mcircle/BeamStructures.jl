function pullback_initialize_beam(grad_output_vec, output_vec, inputs)
    # inputs are node, beam, parameters
    node, beam, parameters = inputs

    # Extract values used in forward pass (or recompute if needed)
    x, y, θ0 = node.x, node.y, node.ϕ
    l, θs, κ = beam.l, beam.θs, beam.κ0
    m_in, fx_in, fy_in = parameters

    # Recompute intermediate values if not captured by AD system
    m_scaled = m_in * normfactor_m(beam)
    fx_scaled = fx_in * normfactor_f(beam)
    fy_scaled = fy_in * normfactor_f(beam)

    # Initialize gradients for inputs
    # Assuming node and beam are structs and we need gradients for their fields
    grad_node_x = zero(typeof(x))
    grad_node_y = zero(typeof(y))
    grad_node_ϕ = zero(typeof(θ0))
    grad_beam_l = zero(typeof(l))
    grad_beam_θs = zero(typeof(θs))
    grad_beam_κ0 = zero(typeof(κ))
    grad_parameters = zero.(parameters) # Gradient for the input parameters vector

    # Assign gradients from grad_output_vec to corresponding output elements' gradients
    grad_m_scaled = grad_output_vec[1]
    grad_θ0_plus_θs = grad_output_vec[2]
    grad_x_div_l = grad_output_vec[3]
    grad_y_div_l = grad_output_vec[4]
    grad_fx_scaled = grad_output_vec[5]
    grad_fy_scaled = grad_output_vec[6]
    grad_κ_times_l = grad_output_vec[7]

    # Backpropagate gradients through the operations

    # Output element 7: κ*l
    grad_κ = grad_κ_times_l * l
    grad_beam_l += grad_κ_times_l * κ
    grad_beam_κ0 += grad_κ # Accumulate to κ0

    # Output element 6: fy_scaled
    grad_fy_in = grad_fy_scaled * normfactor_f(beam)
    # Gradient flows to fy_in (3rd element of parameters)
    grad_parameters[3] += grad_fy_in
    # Gradient flows to normfactor_f(beam)
    grad_normfactor_f_from_fy = grad_fy_scaled * fy_in
    # Propagate grad_normfactor_f_from_fy back through normfactor_f(beam)
    # This requires the pullback of normfactor_f w.r.t. beam parameters.
    # Assuming normfactor_f_pullback(beam)(grad) returns gradient w.r.t. beam fields.
    # Let's represent this conceptually:
    # grad_beam_from_normfactor_f = normfactor_f_pullback(beam)(grad_normfactor_f_from_fy)
    # Accumulate grad_beam_from_normfactor_f to beam gradients.
    # Placeholder: This requires knowledge of normfactor_f implementation.

    # Output element 5: fx_scaled
    grad_fx_in = grad_fx_scaled * normfactor_f(beam)
    # Gradient flows to fx_in (2nd element of parameters)
    grad_parameters[2] += grad_fx_in
    # Gradient flows to normfactor_f(beam)
    grad_normfactor_f_from_fx = grad_fx_scaled * fx_in
    # Accumulate to total grad_normfactor_f
    # grad_normfactor_f_total = grad_normfactor_f_from_fy + grad_normfactor_f_from_fx
    # Propagate grad_normfactor_f_total back to beam (conceptually)

    # Output element 4: y./l
    grad_y = grad_y_div_l / l
    grad_node_y += grad_y
    grad_beam_l += grad_y_div_l * (-y / (l^2))

    # Output element 3: x./l
    grad_x = grad_x_div_l / l
    grad_node_x += grad_x
    grad_beam_l += grad_x_div_l * (-x / (l^2))

    # Output element 2: θ0 + θs
    grad_θ0 = grad_θ0_plus_θs
    grad_node_ϕ += grad_θ0
    grad_θs = grad_θ0_plus_θs
    grad_beam_θs += grad_θs

    # Output element 1: m_scaled
    grad_m_in = grad_m_scaled * normfactor_m(beam)
    # Gradient flows to m_in (1st element of parameters)
    grad_parameters[1] += grad_m_in
    # Gradient flows to normfactor_m(beam)
    grad_normfactor_m = grad_m_scaled * m_in
    # Propagate grad_normfactor_m back through normfactor_m(beam)
    # This requires the pullback of normfactor_m w.r.t. beam parameters.
    # Assuming normfactor_m_pullback(beam)(grad) returns gradient w.r.t. beam fields.
    # Let's represent this conceptually:
    # grad_beam_from_normfactor_m = normfactor_m_pullback(beam)(grad_normfactor_m)
    # Accumulate grad_beam_from_normfactor_m to beam gradients.
    # Placeholder: This requires knowledge of normfactor_m implementation.

    # Accumulate all gradients related to normfactor_f and normfactor_m back to beam parameters
    # This part is conceptual and depends on the implementation of normfactor_m and normfactor_f
    # and their pullbacks. A generic representation:
    # total_grad_normfactor_m_f = grad_normfactor_m + grad_normfactor_f_from_fy + grad_normfactor_f_from_fx
    # grad_beam_from_normfactors = pullback_of_norm_factors_wrt_beam(beam)(total_grad_normfactor_m_f)
    # Add grad_beam_from_normfactors to relevant grad_beam fields.

    # For this example, we can only provide the gradients for the scalar results of normfactor functions.
    # The actual gradient propagation into the `beam` structure from these requires more info.

    # Return gradients for (node, beam, parameters)
    # Assuming node and beam gradients are returned as structures with gradients for their fields
    grad_node = (x=grad_node_x, y=grad_node_y, ϕ=grad_node_ϕ) # Example structure
    grad_beam = (l=grad_beam_l, θs=grad_beam_θs, κ0=grad_beam_κ0) # Example structure
    # Note: The gradients related to normfactor_m and normfactor_f are not fully propagated into grad_beam fields here.

    return grad_node, grad_beam, grad_parameters
end