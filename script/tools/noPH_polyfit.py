import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt


def genLog(x, A, K, nu):
    y = A + ((K - A) / ((1 + 0.0001 * np.exp(-0.65*x))**(1/nu)))
    return y


def expFun(x,a,c):
    return a*np.exp(-0.7*x)-c


def hill(x,h):
    return 281 * x/(x+h)


def main():
    x_fit = np.linspace(0,4,num=40000)
    x = np.linspace(0,14,num=140000)
    pot = np.loadtxt('./potentials/noPH_NB_pot.txt', delimiter='\n')
    pot = np.absolute(pot)
    print(pot[0])

    popt, pcov = opt.curve_fit(genLog, x_fit, pot)
    print(popt)
    y_logistic = genLog(x, *popt)  
    y_log_res = pot - genLog(x_fit, *popt)
    y_log_rmse = (scipy.sum(y_log_res**2)/(y_log_res.size-2))**0.5
    print(f"RMSE of general logarithmic fitting function = {y_log_rmse}")


    popt, pcov = opt.curve_fit(expFun, x_fit, pot)
    print(popt)
    y_exp = expFun(x, *popt)
    y_exp_res = pot - expFun(x_fit, *popt)
    y_exp_rmse = (scipy.sum(y_log_res**2)/(y_log_res.size-2))**0.5
    print(f"RMSE of exponential fitting function = {y_exp_rmse}")


    popt, pcov = opt.curve_fit(hill, x_fit, pot)
    print(popt)
    y_hill = hill(x, *popt)
    y_hill_res = pot - hill(x_fit, *popt)
    y_hill_rmse = (scipy.sum(y_hill_res**2)/(y_hill_res.size-2))**0.5
    print(f"RMSE of Hill fitting function = {y_hill_rmse}")


    fig, ax = plt.subplots()
    ax.set_xlim(0,14)
    ax.set_ylim(0,300)
    ax.grid(b=True, linestyle='--', linewidth=0.5, color='gray')
    ax.plot(x_fit,pot, label="Original potential")
    ax.plot(x,y_logistic, label="fitted logistics curve", lw=1, linestyle='--')
    ax.plot(x,y_exp, label="Fitted exponential curve", lw=1, linestyle='-.')
    ax.plot(x,y_hill, label="Fitted Hill curve", lw=1, linestyle=':')

    ax.legend(frameon=False, loc='lower right')
    plt.show()
      

if __name__ == "__main__":
    main()