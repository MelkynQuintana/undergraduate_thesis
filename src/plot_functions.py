import numpy as np
import matplotlib.pyplot as plt
import corner
from scipy.optimize import curve_fit
from getdist import MCSamples, plots
from matplotlib.patches import Rectangle
import re

'''General plot - 2 Parameters'''

def mc_plot2(sampler, save=False, figname='default'):

    samples = sampler.flatchain

    # Etiquetas para el corner plot
    labels = [r'$H_0$']

    plt.figure(figsize=(10, 6))
    corner.corner(samples, show_titles=True, labels=labels, plot_datapoints=True,
                quantiles=[0.16, 0.5, 0.84])
    
    if save:
        plt.savefig(figname + '.png')
    else:
        plt.show()

'''Plot of linear fit'''

def linearfit_plot(ax, data):
    #Fit data
    params, cov = curve_fit(linear_function, data['Distance'], data['Velocity'],
                             bounds=([-np.inf, 0], [np.inf, 1e-5]))

    m, b = params
    #Model values using the fitted function
    Y_model = linear_function(data['Distance'], m, b)

    #Plot
    ax.errorbar(data['Distance'], data['Velocity'],  fmt='ko', ms=5, ecolor='k', elinewidth=1, capthick=1, capsize=4)
    ax.plot(data['Distance'], Y_model, color='red', label='{:.2f} x + {:.2f}'.format(m,b))
    ax.legend()
    ax.set_xlabel('Distance [Mpc]')
    ax.set_ylabel('Velocity [Km/s]')
    ax.set_title('Fitting')


'''Residual of mcmc results'''
def residuals(data, sampler, model):

    samples = sampler.flatchain

    #Extract mean values from resultos of emcee
    omega_m_model = np.percentile(samples[:,0], [50])
    h0_model = np.percentile(samples[:,1], [50])

    #Model of fit, default FlatModel
    dist_model = model.dist(data['zHD'], h0=h0_model, om_m=omega_m_model)

    X_mc = data['zHD']
    Y_mc = data['Distance']
    Y_e_mc = data['Dist_err']

    residuos_mc = Y_mc - dist_model

    # Crear la figura y los subplots: uno arriba (grafico principal) y otro abajo (residuos)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), gridspec_kw={'height_ratios': [3, 1]})

    # Gráfico principal: datos reales vs. modelo
    ax1.errorbar(X_mc, Y_mc, yerr=Y_e_mc, fmt='ko', ms=5, ecolor='k', elinewidth=1, capthick=1, capsize=4)
    ax1.plot(X_mc, dist_model, color='red')
    ax1.set_ylabel(r'$D_L$')
    ax1.legend()
    ax1.grid(True)

    # Gráfico inferior: residuos (diferencia entre datos reales y modelo)
    ax2.scatter(X_mc, residuos_mc, color='g', s=2, label='Residuos')
    ax2.axhline(0, color='black', linestyle='--')  # Línea en y=0 para referencia
    ax2.set_title('Residuos')
    ax2.set_xlabel('Distance [Mpc]')
    ax2.set_ylabel('Residual')
    ax2.grid(True)

    plt.tight_layout()
    plt.show()


def mu_galmass(mu, galmass):

    mu = np.array(mu)
    galmass = np.array(galmass)

    indices_sort = np.argsort(mu)

    mus = mu[indices_sort]
    mus = mus[:len(mus) - (len(mus) % 5)]

    smass = galmass[indices_sort]
    smass = smass[:len(smass) - (len(smass) % 5)]

    median_mu = np.median(mus.reshape(-1, 5), axis=1)
    median_smass = np.median(smass.reshape(-1, 5), axis=1)

    plt.figure(figsize=(9,6))
    #plt.scatter(matched_galaxies1['MU_SH0ES'], np.log10((1e10)*matched_galaxies1['galaxy_smass']))
    plt.scatter(median_mu, np.log10((1e10) * median_smass))
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$M_{GH}$')
    plt.show()


def h0_corr(mu, galmass, by='median', number=5):

    mu = np.array(mu)
    galmass = np.array(galmass)

    indices_sort = np.argsort(mu)

    mus = mu[indices_sort]
    mus = mus[:len(mus) - (len(mus) % number)]

    smass = galmass[indices_sort]
    smass = smass[:len(smass) - (len(smass) % number)]

    if by == 'median':

        median_mu = np.median(mus.reshape(-1, number), axis=1)
        median_smass = np.median(smass.reshape(-1, number), axis=1)

        plt.figure(figsize=(9,6))
        #plt.scatter(matched_galaxies1['MU_SH0ES'], np.log10((1e10)*matched_galaxies1['galaxy_smass']))
        plt.scatter(median_mu, median_smass)
        plt.xlabel(r'$log(M_{GH})$')
        plt.ylabel(r'$H_0$')
        plt.show()

    elif by == 'mean':

        mean_mu = np.mean(mus.reshape(-1, number), axis=1)
        mean_smass = np.mean(smass.reshape(-1, number), axis=1)

        plt.figure(figsize=(9,6))
        #plt.scatter(matched_galaxies1['MU_SH0ES'], np.log10((1e10)*matched_galaxies1['galaxy_smass']))
        plt.scatter(mean_mu, mean_smass)
        plt.xlabel(r'$log(M_{GH})$')
        plt.ylabel(r'$H_0$')
        plt.show()

'''

def make_triangle_plot(samples, labels, colors, outfile, figsize, dpi, dataset_names):
    g = plots.get_subplot_plotter()
    g.settings.figure_legend_frame = False
    g.settings.legend_fontsize = 14
    g.settings.axes_labelsize = 18  # axis label size
    g.settings.axes_fontsize = 16   # tick label size

    # Auto extract param names
    if isinstance(samples, list):
        param_names = [p.name for p in samples[0].getParamNames().names]
    else:
        param_names = [p.name for p in samples.getParamNames().names]

    g.triangle_plot(
        samples, param_names, filled=True,
        legend_labels=dataset_names if isinstance(samples, list) and len(samples) > 1 else None,
        contour_colors=colors,
        line_args=[{'color': c} for c in colors],
        diag1d_kwargs={'colors': colors}
    )

    g.fig.set_size_inches(*figsize)

    if outfile:
        plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()
'''


def summarize_mc_samples(mc_sample, sample_name):
    """
    Print mean ± standard deviation for each parameter in an MCSamples object.
    
    Parameters
    ----------
    mc_sample : MCSamples
        MCSamples object containing the MCMC samples.
    sample_name : str
        Name of the sample for display
    """
    param_names = mc_sample.getParamNames().names
    for name in param_names:
        mean = mc_sample.mean(name)
        err = mc_sample.std(name)
        print(f"{sample_name} - {name}: {mean:.4f} ± {err:.4f}")

def make_triangle_plot(samples, labels, colors, outfile, figsize, dpi, dataset_names):
    g = plots.get_subplot_plotter()
    g.settings.figure_legend_frame = False
    g.settings.legend_fontsize = 14
    g.settings.axes_labelsize = 18  # axis label size
    g.settings.axes_fontsize = 16   # tick label size

    # Auto extract param names
    if isinstance(samples, list):
        param_names = [p.name for p in samples[0].getParamNames().names]
    else:
        param_names = [p.name for p in samples.getParamNames().names]

    g.triangle_plot(
        samples, param_names, filled=True,  # Sin leyenda
        contour_colors=colors,
        line_args=[{'color': c, 'label': ''} for c in colors],
        diag1d_kwargs={'colors': colors}
    )

    g.fig.set_size_inches(*figsize)

    # Add mean values in the upper-right empty space
    samples_list = samples if isinstance(samples, list) else [samples]
    dataset_names_list = dataset_names if isinstance(dataset_names, list) else [dataset_names]
    colors_list = colors if isinstance(colors, (list, tuple)) else [colors]

    if True:  # Siempre ejecutar

        n_params = len(param_names)
        fig = g.fig
        
        # Try to get existing subplot or find empty space
        try:
            ax_text = fig.add_subplot(g.gridspec[0, n_params-1])
            ax_text.axis('off')
        except:
            ax_text = None
        
        if ax_text is not None:
            # Add text to the subplot centered
            y_start = 0.5  # Centrado verticalmente
            y_step = 0.1
            
            y_pos = y_start + (len(samples_list) * (len(param_names) + 1) * y_step) / 2
            
            for i, (sample, name, color) in enumerate(zip(samples_list, dataset_names_list, colors_list)):
                # Dataset name in bold with LaTeX
                # Separar palabras con regex
                name_spaced = re.sub(r'([a-z])([A-Z])', r'\1 \2', name)
                name_spaced = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1 \2', name_spaced)
                ax_text.text(0.5, y_pos, rf"$\mathbf{{{name_spaced}}}$",
                            transform=ax_text.transAxes,
                            fontsize=16, color=color,
                            verticalalignment='top',
                            horizontalalignment='center')
                y_pos -= y_step
                
                # Parameter means with errors
                for param in param_names:
                    mean = sample.mean(param)
                    err = sample.std(param)
                    
                    # Format parameter name with LaTeX subscripts
                    # Convert param names like "om_m" to "\Omega_m", "H0" to "H_0"
                    if 'om' in param.lower():
                        param_latex = r'\Omega_m'
                    elif 'h0' in param.lower() or param == 'H0':
                        param_latex = r'H_0'
                    else:
                        param_latex = param
                    
                    ax_text.text(0.5, y_pos, 
                                rf"${param_latex} = {mean:.4f} \pm {err:.4f}$", 
                                transform=ax_text.transAxes,
                                fontsize=14, color=color,
                                verticalalignment='top',
                                horizontalalignment='center')
                    y_pos -= y_step
                
                y_pos -= y_step * 0.3  # Extra spacing between datasets
        else:
            # Fallback: add text directly to figure in upper right (centered)
            fig_x = 0.73  # Posición horizontal
            fig_y = 0.65  # Centrado verticalmente sobre el subplot inferior
            y_step = 0.035
            
            y_pos = fig_y + (len(samples_list) * (len(param_names) + 1) * y_step) / 2
            
            for i, (sample, name, color) in enumerate(zip(samples, dataset_names, colors)):
                # Dataset name in bold with LaTeX
                fig.text(fig_x, y_pos, rf"$\mathbf{{{name}}}$", 
                        fontsize=16, color=color,
                        verticalalignment='top',
                        horizontalalignment='center')
                y_pos -= y_step
                
                # Parameter means with errors
                for param in param_names:
                    mean = sample.mean(param)
                    err = sample.std(param)
                    
                    # Format parameter name with LaTeX subscripts
                    if 'om' in param.lower():
                        param_latex = r'\Omega_m'
                    elif 'h0' in param.lower() or param == 'H0':
                        param_latex = r'H_0'
                    else:
                        param_latex = param
                    
                    fig.text(fig_x, y_pos, 
                            rf"${param_latex} = {mean:.4f} \pm {err:.4f}$", 
                            fontsize=14, color=color,
                            verticalalignment='top',
                            horizontalalignment='center')
                    y_pos -= y_step
                
                y_pos -= y_step * 0.3  # Extra spacing between datasets

    g.fig.legends = []

    if outfile:
        plt.savefig(outfile, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()

    return g

'''Linear Function'''
def linear_function(x, m, b):
    return m * x + b    