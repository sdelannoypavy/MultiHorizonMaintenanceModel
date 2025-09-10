
function build_P_linear_daily(nb_state,MTBF)

    p = ((nb_state - 1)/MTBF)/365

    P = [1 - p    0       0;
               p        1 - p   0;
               0        p       1]

end