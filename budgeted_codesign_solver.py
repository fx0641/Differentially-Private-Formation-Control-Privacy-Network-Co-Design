import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigvals
from scipy.stats import norm

class BudgetedCoDesignSolver:
    """Budgeted biconvex co-design solver using alternating convex search."""
    
    def __init__(self, adjacency_matrix, eps_max, err_bound, vartheta=1.0, 
                 b=1.0, delta=0.05, gamma=0.1, epsilon_min=1e-5):
        """Initialize solver with graph and problem parameters."""
        
        self.adj = np.array(adjacency_matrix)
        self.N = self.adj.shape[0]
        self.adj = np.maximum(self.adj, self.adj.T)
        
        self.edges = [(i, j) for i in range(self.N) for j in range(i+1, self.N) 
                      if self.adj[i, j] == 1]
        self.num_edges = len(self.edges)
        
        self.eps_max = np.array(eps_max) if np.isscalar(eps_max) else eps_max
        if self.eps_max.shape == ():
            self.eps_max = np.full(self.N, self.eps_max)
        
        self.err_bound = err_bound
        self.vartheta = vartheta
        self.b = b
        self.delta = delta
        self.gamma = gamma
        self.epsilon_min = epsilon_min
        
        print(f"Initialized BudgetedCoDesignSolver:")
        print(f"  Graph: {self.N} nodes, {self.num_edges} edges")
        print(f"  Parameters: εᵣ={err_bound}, vθ={vartheta}, γ={gamma}")
    
    def laplacian_from_weights(self, weights):
        """Compute Laplacian matrix from edge weights."""
        L = np.zeros((self.N, self.N))
        idx = 0
        for i, j in self.edges:
            w = weights[idx]
            L[i, i] += w
            L[j, j] += w
            L[i, j] -= w
            L[j, i] -= w
            idx += 1
        return L
    
    def get_lambda2(self, weights):
        """Get algebraic connectivity (λ₂)."""
        L = self.laplacian_from_weights(weights)
        eigs = np.linalg.eigvalsh(L)
        eigs = np.sort(eigs)
        return eigs[1]
    
    def sigma_func(self, eps):
        """Compute sigma function."""
        K = norm.isf(self.delta)
        return (self.b / (2 * eps)) * (K + np.sqrt(K**2 + 2 * eps))
    
    def weighted_degree(self, weights):
        """Compute weighted degree per node."""
        deg = np.zeros(self.N)
        idx = 0
        for i, j in self.edges:
            w = weights[idx]
            deg[i] += w
            deg[j] += w
            idx += 1
        return deg
    
    def get_g_function(self, weights):
        """Compute g_i(w) per node."""
        deg = self.weighted_degree(weights)
        g = np.zeros(self.N)
        for i in range(self.N):
            w_sq_sum = np.sum([weights[k]**2 for k, (i2, j2) in enumerate(self.edges) if i2 == i or j2 == i])
            g[i] = w_sq_sum - (deg[i]**2) / self.N
        return g
    
    def get_h_function(self, epsilons):
        """Compute h_i(ε_i) = σ_i²(ε_i)."""
        sigmas = self.sigma_func(epsilons)
        return sigmas**2
    
    def solve_budgeted_acs(self, budget, alpha_lambda2=1.0, max_iterations=30, 
                          tolerance=1e-8, disp=True):
        """Solve budgeted biconvex co-design using alternating convex search."""
        
        if disp:
            print(f"Solving budgeted biconvex co-design using ACS...")
            print(f"  Budget B = {budget} (constraint: Tr(L) ≤ {2*budget})")
            print(f"  Alpha (λ₂ weight) = {alpha_lambda2}")
        
        self.budget = budget
        self.alpha_lambda2 = alpha_lambda2
        
        weights = np.ones(self.num_edges) * 0.5
        epsilons = self.eps_max * 0.5
        y = 0.1
        
        iteration_history = []
        prev_objective = float('inf')
        
        for iteration in range(max_iterations):
            if disp:
                print(f"\n  ACS Iteration {iteration + 1}/{max_iterations}")
                print(f"  {'-'*50}")
            
            weights_prev = weights.copy()
            epsilons_prev = epsilons.copy()
            y_prev = y
            
            if disp:
                print(f"  Step 1: Optimizing weights + connectivity (y) with fixed privacy...")
            weights_new, y_new = self._solve_weights_y_subproblem(weights, epsilons, y)
            
            if disp:
                w_step_change = np.linalg.norm(weights_new - weights)
                y_step_change = abs(y_new - y)
                budget_new = np.trace(self.laplacian_from_weights(weights_new))
                print(f"    → w change: {w_step_change:.2e}, y: {y:.4f}→{y_new:.4f} (Δ={y_step_change:.2e})")
                print(f"    → Budget: {np.trace(self.laplacian_from_weights(weights)):.3f}→{budget_new:.3f}")
                
                wy_result = self._last_weights_y_result
                status_symbol = "✅" if wy_result['success'] else "❌"
                print(f"    → Weights+y solver: {status_symbol} {wy_result['iterations']} iters, {wy_result['function_evals']} f-evals")
                if not wy_result['success']:
                    print(f"      Warning: {wy_result['message']}")
            
            weights, y = weights_new, y_new
            
            if disp:
                print(f"  Step 2: Optimizing privacy parameters with fixed weights + connectivity...")
            epsilons_new = self._solve_epsilons_subproblem(weights, epsilons, y)
            
            if disp:
                eps_step_change = np.linalg.norm(epsilons_new - epsilons)
                avg_eps_old = np.mean(epsilons)
                avg_eps_new = np.mean(epsilons_new)
                print(f"    → ε change: {eps_step_change:.2e}, avg ε: {avg_eps_old:.4f}→{avg_eps_new:.4f}")
                
                eps_result = self._last_epsilons_result
                status_symbol = "✅" if eps_result['success'] else "❌"
                print(f"    → Epsilons solver: {status_symbol} {eps_result['iterations']} iters, {eps_result['function_evals']} f-evals")
                if not eps_result['success']:
                    print(f"      Warning: {eps_result['message']}")
            
            epsilons = epsilons_new
            
            current_objective = self.vartheta * np.sum(epsilons**2) - self.alpha_lambda2 * y
            
            weights_change = np.linalg.norm(weights - weights_prev)
            epsilons_change = np.linalg.norm(epsilons - epsilons_prev)
            y_change = abs(y - y_prev)
            objective_change = abs(current_objective - prev_objective)
            
            iteration_info = {
                'iteration': iteration + 1,
                'objective': current_objective,
                'weights_change': weights_change,
                'epsilons_change': epsilons_change,
                'y_change': y_change,
                'objective_change': objective_change,
                'lambda2_actual': self.get_lambda2(weights),
                'y_auxiliary': y,
                'budget_used': np.trace(self.laplacian_from_weights(weights)),
                'weights_y_iters': self._last_weights_y_result.get('iterations', 0),
                'weights_y_fevals': self._last_weights_y_result.get('function_evals', 0),
                'weights_y_success': self._last_weights_y_result.get('success', False),
                'epsilons_iters': self._last_epsilons_result.get('iterations', 0),
                'epsilons_fevals': self._last_epsilons_result.get('function_evals', 0),
                'epsilons_success': self._last_epsilons_result.get('success', False)
            }
            iteration_history.append(iteration_info)
            
            if disp:
                budget_slack = 2*budget - iteration_info['budget_used']
                
                g = self.get_g_function(weights)
                h = self.get_h_function(epsilons)
                d = 1
                c2 = self.gamma * d / self.N
                s = np.full(self.N, 0.1)
                c1 = (self.N-1) / (self.gamma * self.N**2) * np.sum(s**2)
                left = c1 + c2 * np.sum(g * h)
                right = self.err_bound * y * (2 - self.gamma * y)
                error_slack = right - left
                
                aux_slack = iteration_info['lambda2_actual'] - y
                
                eps_utilization = epsilons / self.eps_max * 100
                avg_eps_util = np.mean(eps_utilization)
                
                print(f"    Objective: {current_objective:.6f} (change: {objective_change:.2e})")
                print(f"    λ₂: {iteration_info['lambda2_actual']:.6f}, y: {y:.6f} (gap: {aux_slack:.2e})")
                print(f"    Budget: {iteration_info['budget_used']:.3f}/{2*budget:.3f} (slack: {budget_slack:.2e})")
                print(f"    Error slack: {error_slack:.2e}")
                print(f"    Privacy util: {avg_eps_util:.1f}% (range: [{np.min(eps_utilization):.1f}%, {np.max(eps_utilization):.1f}%])")
                print(f"    Variable changes: w={weights_change:.2e}, ε={epsilons_change:.2e}, y={y_change:.2e}")
                
                tol = 1e-4
                binding = []
                if budget_slack < tol:
                    binding.append("Budget")
                if error_slack < tol:
                    binding.append("Error")
                if aux_slack < tol:
                    binding.append("Auxiliary")
                
                if binding:
                    print(f"    Binding constraints: {', '.join(binding)}")
                else:
                    print(f"    No binding constraints")
                
                print(f"  {'-'*50}")
            
            if (weights_change < tolerance and 
                epsilons_change < tolerance and 
                y_change < tolerance and
                objective_change < tolerance):
                if disp:
                    print(f"  ✅ ACS converged after {iteration + 1} iterations")
                break
            
            prev_objective = current_objective
        
        else:
            if disp:
                print(f"  ⚠️ ACS reached maximum iterations ({max_iterations})")
        
        final_lambda2 = self.get_lambda2(weights)
        final_budget_used = np.trace(self.laplacian_from_weights(weights))
        
        total_wy_iters = sum(info.get('weights_y_iters', 0) for info in iteration_history)
        total_eps_iters = sum(info.get('epsilons_iters', 0) for info in iteration_history)
        total_wy_fevals = sum(info.get('weights_y_fevals', 0) for info in iteration_history)
        total_eps_fevals = sum(info.get('epsilons_fevals', 0) for info in iteration_history)
        
        if disp:
            print(f"\n  Subproblem solver summary:")
            print(f"    Total ACS iterations: {len(iteration_history)}")
            print(f"    Weights+y subproblems: {total_wy_iters} total SLSQP iterations, {total_wy_fevals} function evaluations")
            print(f"    Epsilons subproblems: {total_eps_iters} total SLSQP iterations, {total_eps_fevals} function evaluations")
            print(f"    Average per ACS iteration: {total_wy_iters/len(iteration_history):.1f} (w+y), {total_eps_iters/len(iteration_history):.1f} (ε)")
        
        return {
            'success': True,
            'message': f'ACS completed after {len(iteration_history)} iterations',
            'weights': weights,
            'epsilons': epsilons,
            'y_auxiliary': y,
            'lambda2_actual': final_lambda2,
            'objective': current_objective,
            'budget': budget,
            'budget_used': final_budget_used,
            'alpha_lambda2': alpha_lambda2,
            'method': 'ACS',
            'iterations': len(iteration_history),
            'iteration_history': iteration_history,
            'converged': len(iteration_history) < max_iterations,
            'subproblem_stats': {
                'total_wy_iters': total_wy_iters,
                'total_eps_iters': total_eps_iters,
                'total_wy_fevals': total_wy_fevals,
                'total_eps_fevals': total_eps_fevals
            }
        }
    
    def _solve_weights_y_subproblem(self, weights_init, epsilons_fixed, y_init):
        """Solve weights + y subproblem with fixed epsilons."""
        
        def objective_wy(x):
            weights = x[:self.num_edges]
            y = x[self.num_edges]
            return -self.alpha_lambda2 * y
        
        def constraint_budget(x):
            weights = x[:self.num_edges]
            L = self.laplacian_from_weights(weights)
            trace_L = np.trace(L)
            return 2*self.budget - trace_L
        
        def constraint_error(x):
            weights = x[:self.num_edges]
            y = x[self.num_edges]
            
            g = self.get_g_function(weights)
            h = self.get_h_function(epsilons_fixed)
            
            d = 1
            c2 = self.gamma * d / self.N
            s = np.full(self.N, 0.1)
            c1 = (self.N-1) / (self.gamma * self.N**2) * np.sum(s**2)
            
            left = c1 + c2 * np.sum(g * h)
            right = self.err_bound * y * (2 - self.gamma * y)
            return right - left
        
        def constraint_auxiliary(x):
            weights = x[:self.num_edges]
            y = x[self.num_edges]
            lambda2 = self.get_lambda2(weights)
            return lambda2 - y
        
        def constraint_weights_nonneg(x):
            weights = x[:self.num_edges]
            return weights
        
        def constraint_y_pos(x):
            y = x[self.num_edges]
            return y - self.epsilon_min
        
        constraints = [
            {'type': 'ineq', 'fun': constraint_budget},
            {'type': 'ineq', 'fun': constraint_error},
            {'type': 'ineq', 'fun': constraint_auxiliary},
            {'type': 'ineq', 'fun': constraint_weights_nonneg},
            {'type': 'ineq', 'fun': constraint_y_pos}
        ]
        
        x0 = np.concatenate([weights_init, [y_init]])
        
        res = minimize(objective_wy, x0, constraints=constraints, 
                      method='SLSQP', options={'disp': False, 'maxiter': 1000})
        
        self._last_weights_y_result = {
            'success': res.success,
            'message': res.message,
            'iterations': res.nit if hasattr(res, 'nit') else 'N/A',
            'function_evals': res.nfev if hasattr(res, 'nfev') else 'N/A',
            'objective_value': res.fun if res.success else 'N/A'
        }
        
        if res.success:
            return res.x[:self.num_edges], res.x[self.num_edges]
        else:
            return weights_init, y_init
    
    def _solve_epsilons_subproblem(self, weights_fixed, epsilons_init, y_fixed):
        """Solve epsilons subproblem with fixed weights and y."""
        
        def objective_eps(epsilons):
            return self.vartheta * np.sum(epsilons**2)
        
        def constraint_error(epsilons):
            g = self.get_g_function(weights_fixed)
            h = self.get_h_function(epsilons)
            
            d = 1
            c2 = self.gamma * d / self.N
            s = np.full(self.N, 0.1)
            c1 = (self.N-1) / (self.gamma * self.N**2) * np.sum(s**2)
            
            left = c1 + c2 * np.sum(g * h)
            right = self.err_bound * y_fixed * (2 - self.gamma * y_fixed)
            return right - left
        
        def constraint_eps_max(epsilons):
            return self.eps_max - epsilons
        
        def constraint_eps_pos(epsilons):
            return epsilons - self.epsilon_min
        
        constraints = [
            {'type': 'ineq', 'fun': constraint_error},
            {'type': 'ineq', 'fun': constraint_eps_max},
            {'type': 'ineq', 'fun': constraint_eps_pos}
        ]
        
        res = minimize(objective_eps, epsilons_init, constraints=constraints,
                      method='SLSQP', options={'disp': False, 'maxiter': 1000})
        
        self._last_epsilons_result = {
            'success': res.success,
            'message': res.message,
            'iterations': res.nit if hasattr(res, 'nit') else 'N/A',
            'function_evals': res.nfev if hasattr(res, 'nfev') else 'N/A',
            'objective_value': res.fun if res.success else 'N/A'
        }
        
        if res.success:
            return res.x
        else:
            return epsilons_init
    
    def analyze_solution(self, result, verbose=True):
        """Analyze optimization solution and constraint slacks."""
        if not result['success']:
            print("❌ Optimization failed, cannot analyze solution")
            return
        
        weights = result['weights']
        epsilons = result['epsilons']
        y_aux = result['y_auxiliary']
        lambda2_actual = result['lambda2_actual']
        budget_used = result['budget_used']
        
        if verbose:
            print(f"\nSolution Analysis")
            print(f"-" * 50)
        
        print(f"Success: {result['success']}")
        print(f"Converged: {result['converged']} in {result['iterations']} iterations")
        print(f"λ₂: {lambda2_actual:.6f}")
        print(f"Auxiliary y: {y_aux:.6f}")
        print(f"Budget used: {budget_used:.6f}/{2*self.budget:.6f} ({budget_used/(2*self.budget)*100:.1f}%)")
        print(f"Objective: {result['objective']:.6f}")
        
        budget_slack = 2*self.budget - budget_used
        
        g = self.get_g_function(weights)
        h = self.get_h_function(epsilons)
        d = 1
        c2 = self.gamma * d / self.N
        s = np.zeros(self.N)
        c1 = (self.N-1) / (self.gamma * self.N**2) * np.sum(s**2)
        left = c1 + c2 * np.sum(g * h)
        right = self.err_bound * y_aux * (2 - self.gamma * y_aux)
        error_slack = right - left
        
        aux_slack = lambda2_actual - y_aux
        
        print(f"\nConstraint slacks:")
        print(f"  Budget: {budget_slack:.6f}")
        print(f"  Error:  {error_slack:.6f}")
        print(f"  Auxiliary: {aux_slack:.6f}")
        
        tol = 1e-3
        binding_constraints = []
        if budget_slack < tol:
            binding_constraints.append("Budget")
        if error_slack < tol:
            binding_constraints.append("Error")
        if aux_slack < tol:
            binding_constraints.append("Auxiliary")
        
        if binding_constraints:
            print(f"  Binding constraints: {', '.join(binding_constraints)}")
        
        eps_utilization = epsilons / self.eps_max * 100
        print(f"\nPrivacy parameter utilization:")
        print(f"  Average: {np.mean(eps_utilization):.1f}%")
        print(f"  Range: [{np.min(eps_utilization):.1f}%, {np.max(eps_utilization):.1f}%]")
        
        return {
            'budget_slack': budget_slack,
            'error_slack': error_slack,
            'aux_slack': aux_slack,
            'binding_constraints': binding_constraints,
            'eps_utilization': eps_utilization
        }