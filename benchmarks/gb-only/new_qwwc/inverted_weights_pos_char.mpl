kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10;
et_hat:=[199855597578490-x_0, -y_0^2*z_0^2-a*y_0^2+a*x_0+x_1, 175506713699382637741901518558-x_1, a*x_1+x_2+(-z_0^2-a)*y_1^2-y_0^2*z_1^2, -b*y_0^2+x_0*z_0^2-b*x_0+y_1^2, d*w_0^3+c*z_0^2-x_0*y_0^2+z_1^2, -147712110368091065533992230151614046615392344-x_2, a*x_2-2*z_1^2*y_1^2+x_3+(-z_0^2-a)*y_2^2-y_0^2*z_2^2, (z_0^2-b)*x_1-b*y_1^2+y_2^2+x_0*z_1^2, d*w_1^3+c*z_1^2-x_0*y_1^2-x_1*y_0^2+z_2^2, f*w_0^3-e*z_0^2+w_1^3-x_0*y_0^2, 112060448118268113054463071321431754045906611580754385116524-x_3, a*x_3-3*z_2^2*y_1^2-3*z_1^2*y_2^2+x_4+(-z_0^2-a)*y_3^2-y_0^2*z_3^2, 2*z_1^2*x_1+(z_0^2-b)*x_2-b*y_2^2+y_3^2+x_0*z_2^2, d*w_2^3+c*z_2^2-x_0*y_2^2-2*x_1*y_1^2-x_2*y_0^2+z_3^2, f*w_1^3-e*z_1^2+w_2^3-x_0*y_1^2-x_1*y_0^2, -28811150117144072850287540203600739227372873503137251992028046798578762646-x_4, a*x_4-4*z_3^2*y_1^2-6*z_2^2*y_2^2-4*z_1^2*y_3^2+x_5+(-z_0^2-a)*y_4^2-y_0^2*z_4^2, 3*z_2^2*x_1+3*z_1^2*x_2+(z_0^2-b)*x_3-b*y_3^2+y_4^2+x_0*z_3^2, d*w_3^3+c*z_3^2-x_0*y_3^2-3*x_1*y_2^2-3*x_2*y_1^2-x_3*y_0^2+z_4^2, f*w_2^3-e*z_2^2+w_3^3-x_0*y_2^2-2*x_1*y_1^2-x_2*y_0^2, -140156802524668290274894967854290075784918950699140673354646660098849987104100137689034540-x_5, a*x_5-5*z_4^2*y_1^2-10*z_3^2*y_2^2-10*z_2^2*y_3^2-5*z_1^2*y_4^2+x_6+(-z_0^2-a)*y_5^2-y_0^2*z_5^2, 4*z_3^2*x_1+6*z_2^2*x_2+4*z_1^2*x_3+(z_0^2-b)*x_4-b*y_4^2+y_5^2+x_0*z_4^2, d*w_4^3+c*z_4^2-x_0*y_4^2-4*x_1*y_3^2-6*x_2*y_2^2-4*x_3*y_1^2-x_4*y_0^2+z_5^2, f*w_3^3-e*z_3^2+w_4^3-x_0*y_3^2-3*x_1*y_2^2-3*x_2*y_1^2-x_3*y_0^2, 549084700982196810537305004873783260026063894259366568641537848560928728708241098157821693601313751612508-x_6, a*x_6-6*z_5^2*y_1^2-15*z_4^2*y_2^2-20*z_3^2*y_3^2-15*z_2^2*y_4^2-6*z_1^2*y_5^2+x_7+(-z_0^2-a)*y_6^2-y_0^2*z_6^2, 5*z_4^2*x_1+10*z_3^2*x_2+10*z_2^2*x_3+5*z_1^2*x_4+(z_0^2-b)*x_5-b*y_5^2+y_6^2+x_0*z_5^2, d*w_5^3+c*z_5^2-x_0*y_5^2-5*x_1*y_4^2-10*x_2*y_3^2-10*x_3*y_2^2-5*x_4*y_1^2-x_5*y_0^2+z_6^2, f*w_4^3-e*z_4^2+w_5^3-x_0*y_4^2-4*x_1*y_3^2-6*x_2*y_2^2-4*x_3*y_1^2-x_4*y_0^2, -1535869454335558310040031757194976981355453852308962614128813375139252073970544967013112537535312890632830578581330254358-x_7, a*x_7-7*z_6^2*y_1^2-21*z_5^2*y_2^2-35*z_4^2*y_3^2-35*z_3^2*y_4^2-21*z_2^2*y_5^2-7*z_1^2*y_6^2+x_8+(-z_0^2-a)*y_7^2-y_0^2*z_7^2, 6*z_5^2*x_1+15*z_4^2*x_2+20*z_3^2*x_3+15*z_2^2*x_4+6*z_1^2*x_5+(z_0^2-b)*x_6-b*y_6^2+y_7^2+x_0*z_6^2, d*w_6^3+c*z_6^2-x_0*y_6^2-6*x_1*y_5^2-15*x_2*y_4^2-20*x_3*y_3^2-15*x_4*y_2^2-6*x_5*y_1^2-x_6*y_0^2+z_7^2, f*w_5^3-e*z_5^2+w_6^3-x_0*y_5^2-5*x_1*y_4^2-10*x_2*y_3^2-10*x_3*y_2^2-5*x_4*y_1^2-x_5*y_0^2, 3788779318895469032209120826981052721256588403982962615252369631057024650046235180800063137969027953582179593872255247615242402309712644-x_8, a*x_8-8*z_7^2*y_1^2-28*z_6^2*y_2^2-56*z_5^2*y_3^2-70*z_4^2*y_4^2-56*z_3^2*y_5^2-28*z_2^2*y_6^2-8*z_1^2*y_7^2+x_9+(-z_0^2-a)*y_8^2-y_0^2*z_8^2, 7*z_6^2*x_1+21*z_5^2*x_2+35*z_4^2*x_3+35*z_3^2*x_4+21*z_2^2*x_5+7*z_1^2*x_6+(z_0^2-b)*x_7-b*y_7^2+y_8^2+x_0*z_7^2, d*w_7^3+c*z_7^2-x_0*y_7^2-7*x_1*y_6^2-21*x_2*y_5^2-35*x_3*y_4^2-35*x_4*y_3^2-21*x_5*y_2^2-7*x_6*y_1^2-x_7*y_0^2+z_8^2, f*w_6^3-e*z_6^2+w_7^3-x_0*y_6^2-6*x_1*y_5^2-15*x_2*y_4^2-20*x_3*y_3^2-15*x_4*y_2^2-6*x_5*y_1^2-x_6*y_0^2, -7745036879377999723512473128067435766521651939720273298734926744112527664574098118242608241146769694913129112722741303938978658673249217102939302836032-x_9, z_aux^6-1];
vars:=[x_9, z_8, y_8, x_8, w_7, z_7, y_7, x_7, w_6, z_6, y_6, x_6, w_5, z_5, y_5, x_5, w_4, z_4, y_4, x_4, w_3, z_3, y_3, x_3, w_2, z_2, y_2, x_2, w_1, z_1, y_1, x_1, w_0, z_0, y_0, x_0, z_aux, w_aux, a, b, c, d, e, f];
new_weights:={w_0 = w_0, w_1 = w_1, w_2 = w_2, w_3 = w_3, w_4 = w_4, w_5 = w_5, w_6 = w_6, w_7 = w_7, w_aux = w_aux, x_0 = x_0^3, x_1 = x_1^3, x_2 = x_2^3, x_3 = x_3^3, x_4 = x_4^3, x_5 = x_5^3, x_6 = x_6^3, x_7 = x_7^3, x_8 = x_8^3, x_9 = x_9^3, y_0 = y_0^2, y_1 = y_1^2, y_2 = y_2^2, y_3 = y_3^2, y_4 = y_4^2, y_5 = y_5^2, y_6 = y_6^2, y_7 = y_7^2, y_8 = y_8^2, z_0 = z_0^2, z_1 = z_1^2, z_2 = z_2^2, z_3 = z_3^2, z_4 = z_4^2, z_5 = z_5^2, z_6 = z_6^2, z_7 = z_7^2, z_8 = z_8^2, z_aux = z_aux};
gb:=CodeTools[Usage](Groebner[Basis](subs(new_weights, et_hat), tdeg(op(vars)), characteristic=11863279),output='all');
# {w_0 = w_0, w_1 = w_1, w_2 = w_2, w_3 = w_3, w_4 = w_4, w_5 = w_5, w_6 = w_6, w_7 = w_7, w_aux = w_aux, x_0 = x_0^3, x_1 = x_1^3, x_2 = x_2^3, x_3 = x_3^3, x_4 = x_4^3, x_5 = x_5^3, x_6 = x_6^3, x_7 = x_7^3, x_8 = x_8^3, x_9 = x_9^3, y_0 = y_0^2, y_1 = y_1^2, y_2 = y_2^2, y_3 = y_3^2, y_4 = y_4^2, y_5 = y_5^2, y_6 = y_6^2, y_7 = y_7^2, y_8 = y_8^2, z_0 = z_0^2, z_1 = z_1^2, z_2 = z_2^2, z_3 = z_3^2, z_4 = z_4^2, z_5 = z_5^2, z_6 = z_6^2, z_7 = z_7^2, z_8 = z_8^2, z_aux = z_aux}
quit;