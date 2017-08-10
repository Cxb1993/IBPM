function Dgxvar = Jac_times( var, varold, gold, params )


epsn = 1e-5;

varpert = varold + epsn * var;

gnew = advance( varpert, params);

Dgxvar = ( gnew - gold ) / epsn; 


Dgxv_n = norm( Dgxvar );

gold_n = norm( gold );


res = norm( Dgxvar - gold ) / norm( gold );


display(['----------------------------------------------------------'])
display(['Dgxv norm = ', num2str(Dgxv_n), ', gold norm = ', num2str( gold_n ),', res = ', num2str(res )])
display(['----------------------------------------------------------'])