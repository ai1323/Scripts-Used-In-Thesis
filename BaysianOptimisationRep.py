import numpy as np
import pandas as pd
import GPyOpt
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt
from GPy.kern import Matern52, RBF

def LoadTheData():
    df = pd.read_csv('/content/data.csv')
    GeometryPARAMETERS = df.iloc[:, :12].values
    performance = df.iloc[:, 12].values
    return GeometryPARAMETERS, performance

def ObjectiveFun(x):
    x = x.reshape(1, -1)
    distances = np.linalg.norm(GeometryPARAMETERS - x, axis=1)
    nearest_index = np.argmin(distances)
    return -performance[nearest_index]

def define_parameter_space(GeometryPARAMETERS):
    bounds = []
    for i in range(12):
        bounds.append({'name': f'param_{i}', 'type': "continuous",
                       'domain': (GeometryPARAMETERS[:, i].min(), GeometryPARAMETERS[:, i].max())})
    return bounds

def setup_models_and_kernels():
    SurrogateModels = {
        'Gaussian Process': 'GP',
        'GP Matern52': 'GP',
        'GP RatQuad': 'GP'
    }

    Kernel = {
        'Gaussian Process': RBF(input_dim=12),
        'GP Matern52': Matern52(input_dim=12),
        'GP RatQuad': RBF(input_dim=12) + RBF(input_dim=12) * RBF(input_dim=12)
    }

    AcqusitionModels = ['EI', 'MPI']

    return SurrogateModels, Kernel, AcqusitionModels

def run_optimization(ObjectiveFun, bounds, SurrogateModels, Kernel, AcqusitionModels, performance):
    results = {}
    for model_name, model in SurrogateModels.items():
        for acq_type in AcqusitionModels:
            key = f"{model_name} - {acq_type}"
            print(f"Optimizing {key}")

            optimizer = BayesianOptimization(
                f=ObjectiveFun,
                domain=bounds,
                model_type='GP',
                acquisition_type=acq_type,
                acquisition_jitter=0.05,
                kernel=Kernel[model_name],
                exact_feval=True,
                maximize=False
            )

            max_iter = 100
            max_time = 300

            for i in range(max_iter):
                if i > 0 and i % 10 == 0:
                    optimizer.model.model.optimize_restarts(num_restarts=5)

                optimizer.run_optimization(max_iter=1, max_time=max_time/max_iter)

                if optimizer.fx_opt < -performance.max():
                    break

            results[key] = {
                'best_params': optimizer.x_opt,
                'best_score': -optimizer.fx_opt,
                'all_performances': -np.array(optimizer.Y_best),
                'all_params': optimizer.X
            }

    return results

def plot_optimization_history(results):
    plt.figure(figsize=(20, 10))
    for key, result in results.items():
        plt.plot(result['all_performances'], label=key)
    plt.xlabel('Iteration')
    plt.ylabel('Best Performance')
    plt.title('Optimization History with Different Surrogate Models and Acquisition Functions')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.tight_layout()
    plt.show()

def plot_parameter_evolution(results):
    plt.figure(figsize=(20, 15))
    for i in range(12):
        plt.subplot(4, 3, i+1)
        for key, result in results.items():
            plt.plot(result['all_params'][:,i], label=key)
        plt.title(f'Parameter: {i}')
        plt.xlabel('Iteration')
        plt.ylabel('Value')
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.show()

def print_best_results(results):
    for key, result in results.items():
        print(f"\nBest results for {key}:")
        print(f"Best score: {result['best_score']}")
        print("Best parameters:")
        for i, param in enumerate(result['best_params']):
            print(f"param_{i}: {param}")
        superformula_params = ', '.join([f"{param:.4f}" for param in result['best_params']])
        print(f"SuperFormula3D({superformula_params})")

def plot_final_performances(results):
    plt.figure(figsize=(15, 8))
    final_performances = [result['best_score'] for result in results.values()]
    plt.bar(results.keys(), final_performances)
    plt.ylabel('Best Performance')
    plt.title('Surrogate Model and Acquisition Function of Final Best')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

def plot_exploration_vs_exploitation(results):
    fig, axes = plt.subplots(3, 3, figsize=(20, 20))
    fig.suptitle('Exploration VS Exploitation through different models and acquisition functions', fontsize=16)

    for idx, (key, result) in enumerate(results.items()):
        row = idx // 3
        col = idx % 3
        ax = axes[row, col]
        scatter = ax.scatter(result['all_params'][:,0], result['all_params'][:,1],
                             c=-result['all_performances'], cmap='viridis')
        ax.set_xlabel('Parameter 0')
        ax.set_ylabel('Parameter 1')
        ax.set_title(key)
        fig.colorbar(scatter, ax=ax)

    plt.tight_layout()
    plt.show()

def main():
    GeometryPARAMETERS, performance = LoadTheData()
    bounds = define_parameter_space(GeometryPARAMETERS)
    SurrogateModels, Kernel, AcqusitionModels = setup_models_and_kernels()
    
    results = run_optimization(ObjectiveFun, bounds, SurrogateModels, Kernel, AcqusitionModels, performance)
    
    plot_optimization_history(results)
    plot_parameter_evolution(results)
    print_best_results(results)
    plot_final_performances(results)
    plot_exploration_vs_exploitation(results)

if __name__ == "__main__":
    main()