
using PyPlot
#function to compute  positive predicted valute PPV (i.e. tp-rate)

#x is a vecotr of zeros and ones (sorted according to some measure, usually the J_ij of DCA)

function compute_PPV(x)
    println("the function return the collect(1:length(x)) and the PPV for x")
    m=length(x)
    tp_check=zeros(m)
    for a=1:m
        tp_check[a]=(x[a]!=0)
    end
    tp_cumsum=cumsum(tp_check)

    tp_rate=zeros(m,2)
    for a=1:m
        tp_rate[a,1]=a
        tp_rate[a,2]=tp_cumsum[a]/a
    end
    return tp_rate[:,1], tp_rate[:,2]
end

function plot_PPV(x; 
                  logscale = false)
    x,y = compute_PPV(x)

    if logscale==true
        plot(x,y)
        xlabel("Num of predictions")
        ylabel("PPV")
        grid(alpha=0.3)
    else
        semilogx(x,y)
        xlabel("Num of predictions")
        ylabel("PPV")
        grid(alpha=0.3)
    end
end

    

#function to compute the ROC curve

function compute_ROC(x)
    println("compute_ROC: the function return the tp_rate and the fp_rate")
    #here are the false positive
    e = 1 - x   #sono i fp

    tpr=cumsum(x)/sum(x)
    fpr=cumsum(e)/sum(e)

    #compute the AUC using trapezoidal rule
    auc = 0.0
    for i in 2:length(fpr)
        dx = fpr[i] - fpr[i-1]
        dy = tpr[i] - tpr[i-1]      
        auc += ( (dx*tpr[i-1]) + (0.5*dx*dy) )
    end

    auc = round(auc,4)

    return tpr, fpr, auc
end

function plot_ROC(x)
    println("the funcion plot the ROC curve and return the AUC (area under curve)")
    tpr,fpr,auc  = compute_ROC(x)

    #plot the ROC curve
    fill_between(fpr,tpr, 0,  alpha=0.3, color="red",label="AUC="*string(auc))
    
    #random case
    a=[0,1]
    b=[0,1]
    plot(a, b, "--")
    grid(alpha=0.3)

    legend()
    return auc
end





