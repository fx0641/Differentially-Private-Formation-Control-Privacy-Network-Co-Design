import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from budgeted_codesign_solver import BudgetedCoDesignSolver
import os
import matplot2tikz
from tikz_utils import generate_all_tikz_files

plt.style.use('ggplot')

def run_budgeted_error_bound_comparison():
    """Run budgeted biconvex optimization with different error bounds."""
    
    print("BUDGETED BICONVEX CO-DESIGN - ERROR BOUND COMPARISON")
    print("-" * 50)
    
    # Graph setup
    
    adj_matrix = np.array([
        [0, 1, 0, 1, 0, 0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 1, 1, 0]
    ])
    
    N = adj_matrix.shape[0]
    edges = [(i, j) for i in range(N) for j in range(i+1, N) if adj_matrix[i, j] == 1]
    
    print(f"Graph: {N} nodes, {len(edges)} edges")
    print(f"Edges: {edges}")
    
    # Problem parameters
    eps_max = np.array([0.4, 0.9, 0.55, 0.35, 0.8, 0.45, 0.7, 0.5, 0.52, 0.58])
    lambda2_min = 0.1
    vartheta = 1.0
    b = 1.0
    delta = 0.05
    gamma = 0.1
    alpha_lambda2 = 1.0
    budget = 6.0
    error_bounds = [8, 16, 64]
    
    print(f"Fixed budget B: {budget} (constraint: Tr(L) ≤ {2*budget})")
    print(f"Error bounds: {error_bounds}")
    print(f"Alpha (λ₂ weight): {alpha_lambda2}")
    
    # Run optimizations
    
    results = {}
    
    for err_bound in error_bounds:
        print(f"\n{'-'*50}")
        print(f"Running optimization with error bound εᵣ = {err_bound}")
        print(f"{'-'*50}")
        
        solver = BudgetedCoDesignSolver(
            adjacency_matrix=adj_matrix,
            eps_max=eps_max,
            err_bound=err_bound,
            vartheta=vartheta,
            b=b,
            delta=delta,
            gamma=gamma
        )
        
        try:
            result = solver.solve_budgeted_acs(
                budget=budget,
                alpha_lambda2=alpha_lambda2,
                max_iterations=30,
                tolerance=1e-8,
                disp=True  # Enable detailed ACS output
            )
            
            if result['success']:
                lambda2 = result['lambda2_actual']
                budget_used = result['budget_used']
                objective = result['objective']
                y_aux = result['y_auxiliary']
                
                # Check constraint tightness
                weights = result['weights']
                epsilons = result['epsilons']
                
                # Budget constraint slack
                budget_slack = 2*budget - budget_used
                
                # Error constraint check
                g = solver.get_g_function(weights)
                h = solver.get_h_function(epsilons)
                d = 1
                c2 = solver.gamma * d / solver.N
                s = np.zeros(solver.N)
                c1 = (solver.N-1) / (solver.gamma * solver.N**2) * np.sum(s**2)
                left = c1 + c2 * np.sum(g * h)
                right = solver.err_bound * y_aux * (2 - solver.gamma * y_aux)
                error_slack = right - left
                
                # Auxiliary constraint slack
                aux_slack = lambda2 - y_aux
                
                print(f"✅ Success! λ₂ = {lambda2:.4f}, Budget used = {budget_used:.2f}/{2*budget:.2f} ({budget_used/(2*budget)*100:.1f}%)")
                print(f"   Objective = {objective:.4f}")
                print(f"   Constraint slacks: Budget={budget_slack:.6f}, Error={error_slack:.6f}, Aux={aux_slack:.6f}")
                
                # Identify binding constraint
                tolerance = 1e-3
                binding_constraints = []
                if budget_slack < tolerance:
                    binding_constraints.append("Budget")
                if error_slack < tolerance:
                    binding_constraints.append("Error")
                if aux_slack < tolerance:
                    binding_constraints.append("Auxiliary")
                
                if binding_constraints:
                    print(f"   Binding constraints: {', '.join(binding_constraints)}")
                else:
                    print(f"   No binding constraints detected")
                
                results[err_bound] = result
                
            else:
                print(f"❌ Failed: {result['message']}")
                results[err_bound] = None
                
        except Exception as e:
            print(f"❌ Error: {e}")
            results[err_bound] = None
    
    # Create figures
    print(f"\nGenerating comparison figures...")
    os.makedirs('figures', exist_ok=True)
    
    valid_results = {er: r for er, r in results.items() if r is not None and r['success']}
    
    if not valid_results:
        print("No successful optimizations to visualize!")
        return
    
    print(f"Creating figures for {len(valid_results)} successful runs...")
    
    create_lambda2_error_figure(valid_results, error_bounds, budget)
    create_privacy_error_figure(valid_results, error_bounds)
    create_edge_weights_error_figure(valid_results, error_bounds)
    create_network_error_visualizations(valid_results, error_bounds, adj_matrix)
    generate_all_tikz_files(valid_results, "$\\varepsilon_r$", "error_bound_comparison", "tikz_budgeted_error")
    
    print(f"\nAnalysis complete!")
    print(f"\nSummary (Fixed Budget B = {budget}):")
    
    temp_solver = BudgetedCoDesignSolver(
        adjacency_matrix=adj_matrix,
        eps_max=eps_max,
        err_bound=1.0,
        vartheta=vartheta,
        b=b,
        delta=delta,
        gamma=gamma
    )
    
    for err_bound, result in sorted(valid_results.items()):
        if result:
            lambda2 = result['lambda2_actual']
            budget_used = result['budget_used']
            utilization = budget_used / (2*budget) * 100
            avg_epsilon = np.mean(result['epsilons'])
            objective = result['objective']
            
            g = temp_solver.get_g_function(result['weights'])
            h = temp_solver.get_h_function(result['epsilons'])
            d = 1
            c2 = temp_solver.gamma * d / temp_solver.N
            s = np.full(temp_solver.N, 0.1)
            c1 = (temp_solver.N-1) / (temp_solver.gamma * temp_solver.N**2) * np.sum(s**2)
            steady_state_error = (c1 + c2 * np.sum(g * h)) / (temp_solver.N * lambda2 * (2 - temp_solver.gamma * lambda2))
            
            print(f"εᵣ={err_bound:2d}: λ₂={lambda2:.4f}, Used={budget_used:5.2f}/{2*budget:5.2f} ({utilization:5.1f}%), Avg ε={avg_epsilon:.3f}, Obj={objective:.3f}, SS Error={steady_state_error:.4f}")
    
    print(f"\nOutput files:")
    print(f"  - tikz_budgeted_error/error_bound_comparison_*.tex")

def create_lambda2_error_figure(valid_results, error_bounds, budget):
    """Plot λ₂ vs error bound."""
    
    error_bounds_list = []
    lambda2_list = []
    budget_used_list = []
    
    for err_bound in sorted(valid_results.keys()):
        result = valid_results[err_bound]
        error_bounds_list.append(err_bound)
        lambda2_list.append(result['lambda2_actual'])
        budget_used_list.append(result['budget_used'])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # λ₂ vs Error Bound - use regular plot if all values are the same
    if len(set(lambda2_list)) == 1:
        # All lambda2 values are the same, use regular plot
        ax1.plot(error_bounds_list, lambda2_list, 'o-', linewidth=2, markersize=8, color='#2E86AB')
        ax1.set_xscale('log')
    else:
        ax1.semilogx(error_bounds_list, lambda2_list, 'o-', linewidth=2, markersize=8, color='#2E86AB')
    
    ax1.set_xlabel('Error Bound (εᵣ)', fontsize=12)
    ax1.set_ylabel('Algebraic Connectivity (λ₂)', fontsize=12)
    ax1.set_title(f'λ₂ vs Error Bound (Budget B = {budget})', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for x, y in zip(error_bounds_list, lambda2_list):
        ax1.text(x, y + 0.01, f'{y:.4f}', ha='center', va='bottom', fontsize=10)
    
    # Budget utilization vs error bound
    utilization = [used/(2*budget)*100 for used in budget_used_list]
    
    # Use regular plot with log scale to avoid rendering issues
    ax2.plot(error_bounds_list, utilization, 's-', linewidth=2, markersize=8, color='#A23B72')
    ax2.set_xscale('log')
    ax2.set_xlabel('Error Bound (εᵣ)', fontsize=12)
    ax2.set_ylabel('Budget Utilization (%)', fontsize=12)
    ax2.set_title(f'Budget Utilization vs Error Bound', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add utilization values
    for x, y in zip(error_bounds_list, utilization):
        ax2.text(x, y + 1, f'{y:.1f}%', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.show()

def create_privacy_error_figure(valid_results, error_bounds):
    """Plot privacy parameters vs error bounds."""
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = ['#2E86AB', '#A23B72', '#F18F01']
    node_indices = np.arange(10)
    bar_width = 0.25
    
    for i, err_bound in enumerate(sorted(valid_results.keys())):
        result = valid_results[err_bound]
        epsilons = result['epsilons']
        
        offset = (i - len(valid_results)/2 + 0.5) * bar_width
        bars = ax.bar(node_indices + offset, epsilons, bar_width, 
                     label=f'e_r = {err_bound}', color=colors[i % len(colors)], alpha=0.8)
    
    ax.set_xlabel('Node Index', fontsize=12)
    ax.set_ylabel('Privacy Parameter (e)', fontsize=12)
    ax.set_title('Privacy Parameters vs Error Bound\n(Budgeted Biconvex Co-design)', fontsize=14, fontweight='bold')
    ax.set_xticks(node_indices)
    ax.set_xticklabels([f'Node {i+1}' for i in range(10)])
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.show()

def create_edge_weights_error_figure(valid_results, error_bounds):
    """Plot edge weights vs error bounds."""
    
    n_edges = len(next(iter(valid_results.values()))['weights'])
    
    fig, ax = plt.subplots(figsize=(16, 6))
    
    colors = ['#2E86AB', '#A23B72', '#F18F01']
    edge_indices = np.arange(n_edges)
    bar_width = 0.25
    
    for i, err_bound in enumerate(sorted(valid_results.keys())):
        result = valid_results[err_bound]
        weights = result['weights']
        
        offset = (i - len(valid_results)/2 + 0.5) * bar_width
        ax.bar(edge_indices + offset, weights, bar_width, 
               label=f'εᵣ = {err_bound}', color=colors[i % len(colors)], alpha=0.8)
    
    ax.set_xlabel('Edge Index', fontsize=12)
    ax.set_ylabel('Edge Weight', fontsize=12)
    ax.set_title('Edge Weights vs Error Bound\n(Budgeted Biconvex Co-design)', fontsize=14, fontweight='bold')
    ax.set_xticks(edge_indices)
    ax.set_xticklabels([f'E{i+1}' for i in range(n_edges)], rotation=45)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.show()

def create_network_error_visualizations(valid_results, error_bounds, adj_matrix):
    """Visualize networks for different error bounds."""
    
    n_plots = len(valid_results)
    fig, axes = plt.subplots(1, n_plots, figsize=(6*n_plots, 6))
    
    if n_plots == 1:
        axes = [axes]
    
    # Create NetworkX graph structure (same as privacy_parameter_comparison.py)
    N = adj_matrix.shape[0]
    edges = [(i, j) for i in range(N) for j in range(i+1, N) if adj_matrix[i, j] == 1]
    
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)
    pos = nx.circular_layout(G)  # Same as original: circular layout
    
    # Calculate global vmin/vmax for consistent coloring across all plots
    all_epsilons = []
    for result in valid_results.values():
        all_epsilons.extend(result['epsilons'])
    global_vmin = min(all_epsilons)
    global_vmax = max(all_epsilons)
    
    # Plot each result
    for idx, err_bound in enumerate(sorted(valid_results.keys())):
        result = valid_results[err_bound]
        ax = axes[idx]
        
        weights = result['weights']
        epsilons = result['epsilons']
        lambda2 = result['lambda2_actual']
        objective = result['objective']
        
        # Draw nodes with color mapped to epsilon (same as original)
        cmap = plt.get_cmap('viridis')
        nx.draw_networkx_nodes(G, pos, node_color=epsilons.tolist(), 
                              node_size=800, cmap=cmap, vmin=global_vmin, vmax=global_vmax, ax=ax)
        
        # Draw edges with width mapped to weights (same scaling as original)
        edge_widths = [max(0.5, w*3) for w in weights]  # Same scaling as original
        for i, (u, v) in enumerate(edges):
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], 
                                 width=edge_widths[i], ax=ax, alpha=0.7)
        
        # Draw node labels (privacy parameters) - same style as original
        node_labels = {i: f'{epsilons[i]:.2f}' for i in range(10)}
        nx.draw_networkx_labels(G, pos, labels=node_labels, 
                               font_color='white', font_size=10, ax=ax)
        
        # Draw edge labels (weights) - same style as original
        edge_labels = {(i, j): f'{weights[k]:.1f}' for k, (i, j) in enumerate(edges)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, 
                                   font_color='red', font_size=8, ax=ax)
        
        # Set title and remove axis (same style as original, but for error bound)
        ax.set_title(f'$\\varepsilon_r = {err_bound}$\nObjective = {objective:.4f}', 
                    fontsize=12, fontweight='bold')
        ax.axis('off')
        
        # Add colorbar for the first subplot (same as original)
        if idx == 0:
            from matplotlib import colors
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=colors.Normalize(vmin=global_vmin, vmax=global_vmax))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.7)
            cbar.set_label('Privacy Parameter (ε)', fontsize=10)
    
    plt.suptitle('Network Graphs with Privacy Parameters and Edge Weights\n(Budgeted Biconvex Co-design)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_budgeted_error_bound_comparison()