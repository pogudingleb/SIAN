infolevel[Groebner]:=10;
Et_hat := [2191010808179183975-w_0, c_0*q*w_0*y_0-c_0*w_0*x_0*y_0+b*w_0+w_1, 424845563503305741825-z_0, -c_0*q*w_0*y_0+h*z_0+z_1, 105973668278158569615-beta_0, beta_1, 290231919917271967234-c_0, c_1, -w_1+4576762280351651251210437785546766549869514815025845657628998418653608771242925, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c_0+y_0*c_1*(q-x_0)*w_0+b*w_1+w_2, beta_0*v_0*x_0+d*x_0-lm+x_1, -beta_0*v_0*x_0+a*y_0+y_1, 92737915633146731106512582122603697213400333703326646233014715158484701288150-z_1, ((-w_0*y_1-w_1*y_0)*c_0-y_0*c_1*w_0)*q+h*z_1+z_2, -w_2+9560314760956026997523174387845059819902017764853643151721356085699311487277515431326363266425589070391700219367502475204932306475357139275, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(q*w_2-2*w_1*x_1-w_2*x_0)*y_0-2*y_1*w_1*(x_0-q))*c_0+((-2*c_1*x_1+c_2*q-c_2*x_0)*y_0-2*y_1*c_1*(x_0-q))*w_0-2*c_1*w_1*(x_0-q)*y_0+b*w_2+w_3, c_2, beta_0*x_0*v_1+v_0*x_0*beta_1+(beta_0*v_0+d)*x_1+x_2, (-beta_0*v_1-beta_1*v_0)*x_0-v_0*beta_0*x_1+a*y_1+y_2, -k*y_0+u*v_0+v_1, 193718530572172524116698724067263174200745305121194703124226011907156447333945690256604182412116859325430703844213184681819826456732314800-z_2, ((-w_0*y_2-2*w_1*y_1-w_2*y_0)*c_0+(-2*c_1*y_1-c_2*y_0)*w_0-2*c_1*w_1*y_0)*q+h*z_2+z_3, -w_3+19970366108141210123349792360301842449761177204471298903151265526227453957301301696800405399730419921696829656544747899290157189922282786897186167646005585264024859279574066836310314039053628369779575, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(-3*w_1*y_2-3*w_2*y_1)*x_0+(3*w_1*y_2+3*w_2*y_1)*q-6*x_1*y_1*w_1)*c_0+((-3*c_1*x_2-3*c_2*x_1+c_3*q-c_3*x_0)*y_0+(-3*c_1*y_2-3*c_2*y_1)*x_0+(3*c_1*y_2+3*c_2*y_1)*q-6*x_1*y_1*c_1)*w_0+((-3*c_1*w_2-3*c_2*w_1)*x_0+(3*c_1*w_2+3*c_2*w_1)*q-6*x_1*c_1*w_1)*y_0-6*c_1*w_1*x_0*y_1+6*c_1*q*w_1*y_1+b*w_3+w_4, c_3, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta_0+(2*beta_1*x_1+beta_2*x_0)*v_0+2*beta_1*v_1*x_0+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta_0+(-2*beta_1*x_1-beta_2*x_0)*v_0-2*beta_1*v_1*x_0+a*y_2+y_3, beta_2, -k*y_1+u*v_1+v_2, 404655084501691642598287607053420661582227482635223783405190834388578369392942690976817858548934607862845207910801527257194325332627265683077613084071098374794574912348500318033498507674578641780850-z_3, ((-w_0*y_3-3*w_1*y_2-3*w_2*y_1-w_3*y_0)*c_0+(-3*c_1*y_2-3*c_2*y_1-c_3*y_0)*w_0+(-3*c_1*w_2-3*c_2*w_1)*y_0-6*c_1*w_1*y_1)*q+h*z_3+z_4, -w_4+41715731381768202156949040644489514799263577252827873802848004904142493094185062608256858986668790366817064328457309757106447197377382219062448003051275331891563349253537086074400668158276294817357133895254011174428072098816009582886325340425353103018927650475, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(-4*w_1*y_3-6*w_2*y_2-4*w_3*y_1)*x_0+(4*w_1*y_3+6*w_2*y_2+4*w_3*y_1)*q+(-12*x_1*y_2-12*x_2*y_1)*w_1-12*w_2*x_1*y_1)*c_0+((-4*c_1*x_3-6*c_2*x_2-4*c_3*x_1+c_4*q-c_4*x_0)*y_0+(-4*c_1*y_3-6*c_2*y_2-4*c_3*y_1)*x_0+(4*c_1*y_3+6*c_2*y_2+4*c_3*y_1)*q+(-12*x_1*y_2-12*x_2*y_1)*c_1-12*x_1*y_1*c_2)*w_0+((-4*c_1*w_3-6*c_2*w_2-4*c_3*w_1)*x_0+(4*c_1*w_3+6*c_2*w_2+4*c_3*w_1)*q+(-12*w_1*x_2-12*w_2*x_1)*c_1-12*x_1*c_2*w_1)*y_0+((-12*w_1*y_2-12*w_2*y_1)*c_1-12*y_1*c_2*w_1)*x_0+((12*w_1*y_2+12*w_2*y_1)*c_1+12*y_1*c_2*w_1)*q-24*c_1*w_1*x_1*y_1+b*w_4+w_5, c_4, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta_0+(3*beta_1*x_2+3*beta_2*x_1+beta_3*x_0)*v_0+(3*beta_1*v_2+3*beta_2*v_1)*x_0+6*beta_1*v_1*x_1+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta_0+(-3*beta_1*x_2-3*beta_2*x_1-beta_3*x_0)*v_0+(-3*beta_1*v_2-3*beta_2*v_1)*x_0-6*beta_1*v_1*x_1+a*y_3+y_4, beta_3, -k*y_2+u*v_2+v_3, 845276582108212199688726213843667502621010296746811480191733820639505352710417007170821730560167477352445132606068197064853673313261879114170247913681660275949880274167263873538993544201676831110130305584547888719281784538949615123039955862427994518950695200-z_4, ((-w_0*y_4-4*w_1*y_3-6*w_2*y_2-4*w_3*y_1-w_4*y_0)*c_0+(-4*c_1*y_3-6*c_2*y_2-4*c_3*y_1-c_4*y_0)*w_0+(-4*c_1*w_3-6*c_2*w_2-4*c_3*w_1)*y_0+(-12*w_1*y_2-12*w_2*y_1)*c_1-12*y_1*c_2*w_1)*q+h*z_4+z_5, -w_5+87139225955723548222161376501541473453661358007082321828606203321057011029283148705129445576038931623923796244534293915237593522203506781203746641610551348413992232847792927716931175444064430721421194082506546865139560416017278858495304971929428281720675225389586943720905009833117489653315079757352956916602474174876175, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(-5*w_1*y_4-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1)*q+(-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-30*w_2*x_2-20*w_3*x_1)*y_1-30*x_1*y_2*w_2)*c_0+((-5*c_1*x_4-10*c_2*x_3-10*c_3*x_2-5*c_4*x_1+c_5*q-c_5*x_0)*y_0+(-5*c_1*y_4-10*c_2*y_3-10*c_3*y_2-5*c_4*y_1)*x_0+(5*c_1*y_4+10*c_2*y_3+10*c_3*y_2+5*c_4*y_1)*q+(-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*c_1+(-30*c_2*x_2-20*c_3*x_1)*y_1-30*x_1*y_2*c_2)*w_0+((-5*c_1*w_4-10*c_2*w_3-10*c_3*w_2-5*c_4*w_1)*x_0+(5*c_1*w_4+10*c_2*w_3+10*c_3*w_2+5*c_4*w_1)*q+(-20*w_1*x_3-30*w_2*x_2-20*w_3*x_1)*c_1+(-30*c_2*x_2-20*c_3*x_1)*w_1-30*w_2*x_1*c_2)*y_0+((-20*w_1*y_3-30*w_2*y_2-20*w_3*y_1)*c_1+(-30*c_2*y_2-20*c_3*y_1)*w_1-30*w_2*y_1*c_2)*x_0+((20*w_1*y_3+30*w_2*y_2+20*w_3*y_1)*c_1+(30*c_2*y_2+20*c_3*y_1)*w_1+30*w_2*y_1*c_2)*q+((-60*x_1*y_2-60*x_2*y_1)*w_1-60*w_2*x_1*y_1)*c_1-60*c_2*w_1*x_1*y_1+b*w_5+w_6, c_5, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta_0+(4*beta_1*x_3+6*beta_2*x_2+4*beta_3*x_1+beta_4*x_0)*v_0+(4*beta_1*v_3+6*beta_2*v_2+4*beta_3*v_1)*x_0+(12*v_1*x_2+12*v_2*x_1)*beta_1+12*beta_2*v_1*x_1+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta_0+(-4*beta_1*x_3-6*beta_2*x_2-4*beta_3*x_1-beta_4*x_0)*v_0+(-4*beta_1*v_3-6*beta_2*v_2-4*beta_3*v_1)*x_0+(-12*v_1*x_2-12*v_2*x_1)*beta_1-12*beta_2*v_1*x_1+a*y_4+y_5, beta_4, -k*y_3+u*v_3+v_4, 1765682744701937129521345737288847917894673872114696628585230626627470154862202021030776232352271262409206418184170700980702791003897854849086576432241116231031158238528908202501196497374621529482839060877454577823065278835374231700995691338749136145582397842121359451589292932705832056118126341476246209399586302317650-z_5, ((-w_0*y_5-5*w_1*y_4-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1-w_5*y_0)*c_0+(-5*c_1*y_4-10*c_2*y_3-10*c_3*y_2-5*c_4*y_1-c_5*y_0)*w_0+(-5*c_1*w_4-10*c_2*w_3-10*c_3*w_2-5*c_4*w_1)*y_0+(-20*w_1*y_3-30*w_2*y_2-20*w_3*y_1)*c_1+(-30*c_2*y_2-20*c_3*y_1)*w_1-30*w_2*y_1*c_2)*q+h*z_5+z_6, -w_6+182023530420019454235242955535250941680125631096576461946289814872862921039150394607778231983943507080848626240204826648979664637509282378606341842343673565639335726806704678402143341188694224637281440285234300883135498627758680549897361245393235136326808484385342714104248725079012016246943857449610524506829249830226652432850272990938812265162973259190495857229333921879375127775, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(-6*w_1*y_5-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*w_1*y_5+15*w_2*y_4+20*w_3*y_3+15*w_4*y_2+6*w_5*y_1)*q+(-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(-60*w_2*y_3-60*w_3*y_2)*x_1-90*w_2*x_2*y_2)*c_0+((-6*c_1*x_5-15*c_2*x_4-20*c_3*x_3-15*c_4*x_2-6*c_5*x_1+c_6*q-c_6*x_0)*y_0+(-6*c_1*y_5-15*c_2*y_4-20*c_3*y_3-15*c_4*y_2-6*c_5*y_1)*x_0+(6*c_1*y_5+15*c_2*y_4+20*c_3*y_3+15*c_4*y_2+6*c_5*y_1)*q+(-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*c_1+(-60*c_2*x_3-60*c_3*x_2-30*c_4*x_1)*y_1+(-60*c_2*y_3-60*c_3*y_2)*x_1-90*x_2*y_2*c_2)*w_0+((-6*c_1*w_5-15*c_2*w_4-20*c_3*w_3-15*c_4*w_2-6*c_5*w_1)*x_0+(6*c_1*w_5+15*c_2*w_4+20*c_3*w_3+15*c_4*w_2+6*c_5*w_1)*q+(-30*w_1*x_4-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*c_1+(-60*c_2*x_3-60*c_3*x_2-30*c_4*x_1)*w_1+(-60*c_2*w_3-60*c_3*w_2)*x_1-90*w_2*x_2*c_2)*y_0+((-30*w_1*y_4-60*w_2*y_3-60*w_3*y_2-30*w_4*y_1)*c_1+(-60*c_2*y_3-60*c_3*y_2-30*c_4*y_1)*w_1+(-60*c_2*w_3-60*c_3*w_2)*y_1-90*w_2*y_2*c_2)*x_0+((30*w_1*y_4+60*w_2*y_3+60*w_3*y_2+30*w_4*y_1)*c_1+(60*c_2*y_3+60*c_3*y_2+30*c_4*y_1)*w_1+(60*c_2*w_3+60*c_3*w_2)*y_1+90*w_2*y_2*c_2)*q+((-120*x_1*y_3-180*x_2*y_2-120*x_3*y_1)*w_1+(-180*w_2*x_2-120*w_3*x_1)*y_1-180*x_1*y_2*w_2)*c_1+((-180*c_2*x_2-120*c_3*x_1)*y_1-180*x_1*y_2*c_2)*w_1-180*c_2*w_2*x_1*y_1+b*w_6+w_7, c_6, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta_0+(5*beta_1*x_4+10*beta_2*x_3+10*beta_3*x_2+5*beta_4*x_1+beta_5*x_0)*v_0+(5*beta_1*v_4+10*beta_2*v_3+10*beta_3*v_2+5*beta_4*v_1)*x_0+(20*v_1*x_3+30*v_2*x_2+20*v_3*x_1)*beta_1+(30*beta_2*x_2+20*beta_3*x_1)*v_1+30*beta_2*v_2*x_1+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta_0+(-5*beta_1*x_4-10*beta_2*x_3-10*beta_3*x_2-5*beta_4*x_1-beta_5*x_0)*v_0+(-5*beta_1*v_4-10*beta_2*v_3-10*beta_3*v_2-5*beta_4*v_1)*x_0+(-20*v_1*x_3-30*v_2*x_2-20*v_3*x_1)*beta_1+(-30*beta_2*x_2-20*beta_3*x_1)*v_1-30*beta_2*v_2*x_1+a*y_5+y_6, beta_5, -k*y_4+u*v_4+v_5, -w_7+380225613243372033902326833293793384609755437660034696518048448810291809720013044584631492259226102706194564293910231600269067985514361263776662272646499026686221217926044162151469979907878968945868040454268057838567915871074230414172513726399795248865682335490483683762544474311233514302269983586336213293245739689843649272291985696347021219592623657367053244027869198922173860161696663534342531869041103362396423692866792864728712156844075, 3688302291733248072456014844465909234029440425616349540401011452478870314516348409495576209570545583552941940143066289067381556663488529690418365018007942959473823568755272879331291513859886342556270357613886178303839096465615375352613413824838998366871557134563696861784894487362772380714647858703951155694319115702072395103838848663158488748040979020142080431705318844599924800-z_6, -beta_1, -beta_2, -beta_3, -beta_4, -beta_5, -c_1, -c_2, -c_3, -c_4, -c_5, -c_6, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, c_6, beta_5, z_5, y_5, x_5, w_5, v_5, c_5, beta_4, z_4, y_4, x_4, w_4, v_4, c_4, beta_3, z_3, y_3, x_3, w_3, v_3, c_3, beta_2, z_2, y_2, x_2, w_2, v_2, c_2, beta_1, z_1, y_1, x_1, w_1, v_1, c_1, beta_0, z_0, y_0, x_0, w_0, v_0, c_0, z_aux, w_aux, a, b, d, h, k, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [beta, c] 658
# 266, [1 = 8, q*w_3*y_1 = 1, b*w_2 = 1, a*y_0 = 1, -5*w_4*x_0*y_1 = 1, -w_0*x_0*y_1 = 1, -z_6 = 1, -4*w_0*x_3*y_1 = 1, v_2*x_2 = 1, -5*w_1*x_0*y_4 = 1, -q*w_0*y_3 = 1, -5*v_1*x_4 = 1, u*v_2 = 1, -2*y_1*w_1*q = 1, -6*x_1*y_1*w_1 = 1, -20*w_3*x_3*y_0 = 1, -x_3*v_0 = 1, -w_0*x_0*y_2 = 1, -q*w_0*y_2 = 1, -k*y_2 = 1, -10*v_2*x_3 = 1, -3*q*w_1*y_2 = 1, -15*w_4*x_0*y_2 = 1, -60*w_1*x_2*y_3 = 1, q*w_3*y_3 = 1, v_3*x_0 = 1, -3*w_1*x_2*y_0 = 1, w_1 = 1, h*z_2 = 1, w_1*y_0*q = 1, u*v_3 = 1, x_6 = 1, d*x_3 = 1, -164721582855644909701359059267008809550018413877809439532012416060209601998657816842945745567897204576735231792678462717632891988364645853664005984739977361895622179151588377402130131036243839382931271497947324584784195495108097455350355250611728672076597977632656 = 1, -w_0*x_4*y_0 = 1, v_3*x_2 = 1, z_6 = 1, -2*w_0*x_1*y_1 = 1, -v_0*x_5 = 1, -60*x_1*y_2*w_3 = 1, -q*w_0*y_0 = 1, v_5 = 1, -q*w_4*y_0 = 1, -60*w_1*x_3*y_2 = 1, y_2*q*w_2 = 1, -w_4*x_0*y_0 = 1, v_2*x_0 = 1, -10*w_3*x_2*y_0 = 1, x_2 = 1, -3*x_2*v_1 = 1, b*w_6 = 1, -4*w_1*x_0*y_3 = 1, z_1 = 1, -20*w_3*x_1*y_1 = 1, u*v_4 = 1, v_4 = 1, -20*w_1*x_1*y_3 = 1, a*y_2 = 1, -3*w_2*y_1*x_0 = 1, v_0*x_0 = 1, -30*w_4*x_1*y_1 = 1, -4*q*w_3*y_1 = 1, -z_2 = 1, q*w_2*y_3 = 1, h*z_5 = 1, -z_5 = 1, y_6 = 1, -60*w_3*x_2*y_1 = 1, -3811631523823849533748331710598418339419473585110868329611049199088727487210518358698403600650604163616865365446917174642824731965169225010856 = 1, -20*w_3*x_0*y_3 = 1, q*w_4*y_0 = 1, -k*y_0 = 1, -q*w_0*y_5 = 1, -815402048292665350361605655594682735739551400388698172957514912181357529452166069241702849357162652363127320586148655376530353341138636285780490937666199100656838591589818068069583942858123853509856470728047200160404503620011204187998022455209975923881995078157740323056620551393441867329220664704405892817577299688714602372 = 1, -w_0*x_6*y_0 = 1, v_2*x_3 = 1, -5*v_4*x_1 = 1, x_3 = 1, -q*w_3*y_0 = 1, q*w_3*y_2 = 1, -15*w_4*x_2*y_0 = 1, -w_5 = 1, -3*w_2*x_1*y_0 = 1, z_4 = 1, -v_3*x_0 = 1, -10*w_3*x_0*y_2 = 1, q*w_0*y_3 = 1, q*w_0*y_6 = 1, -7118526460094158566237096960490561709091070092828086062441195468912695928415522698155996215931724178282067881868879638064579962114818990405168966475255010638754980655611998225896646171876010907661249599240358565428108238285506433877468031680719737163294282252626338152191118292372037938699102061638984825319150824776777527388396550000914304469622110278259228804455923450618434869065040 = 1, -3*v_2*x_1 = 1, -10*w_0*x_3*y_2 = 1, b*w_1 = 1, v_0*x_5 = 1, -w_1*y_0*x_0 = 1, y_4 = 1, -4*w_3*x_1*y_0 = 1, -w_0*x_0*y_6 = 1, -18868275170634996523136416062498983968867351914014193414471986630184526794799440538880069229480783183353398157120903083393105059065818070513252870661353206952132444109638512536845548742781744068390199300 = 1, -35238011654326437371059416444119838700537121247060979719137014113665350279779162354060118065937189034372584167858558515279505959401029616894137584301857373169751419936952549627317562604131276306642254481919809910306853648744994676069811331658693417913129106079246372861100360171488943149337250394854532631205302046461561062085283811469131525517147743176523043922730426461827142223497192841862388006787187122123531501245381782281308422794350023432 = 1, -4*w_0*x_1*y_3 = 1, -5*q*w_1*y_4 = 1, y_1 = 1, v_3*x_1 = 1, v_2*x_1 = 1, -w_3*x_0*y_0 = 1, -z_1 = 1, a*y_5 = 1, -6*w_5*x_0*y_1 = 1, -6*w_0*x_2*y_2 = 1, -x_0*v_1 = 1, -w_0*x_5*y_0 = 1, v_5*x_0 = 1, b*w_3 = 1, -w_0*x_1*y_0 = 1, x_0*v_1 = 1, -6*w_1*x_0*y_5 = 1, q*w_4*y_2 = 1, h*z_3 = 1, x_4*v_0 = 1, u*v_1 = 1, w_5 = 1, -6*w_2*x_2*y_0 = 1, -10*w_2*x_3*y_0 = 1, q*w_1*y_2 = 1, a*y_1 = 1, -10*q*w_2*y_3 = 1, -30*w_1*x_1*y_4 = 1, v_4*x_0 = 1, -3*w_2*y_1*q = 1, -10*v_3*x_2 = 1, d*x_5 = 1, -4*w_3*x_0*y_1 = 1, v_1*x_4 = 1, -5*w_0*x_4*y_1 = 1, y_5 = 1, -w_7 = 1, -6*w_5*x_1*y_0 = 1, z_5 = 1, b*w_5 = 1, v_1 = 1, h*z_1 = 1, z_aux = 1, -15*w_2*x_4*y_0 = 1, a*y_4 = 1, -q*w_5*y_0 = 1, x_1 = 1, -4*v_1*x_3 = 1, -w_2 = 1, z_2 = 1, w_7 = 1, -60*w_2*x_3*y_1 = 1, -w_4 = 1, -10*w_2*x_0*y_3 = 1, d*x_4 = 1, -w_6*x_0*y_0 = 1, -15*w_0*x_4*y_2 = 1, -90*w_2*x_2*y_2 = 1, b*w_4 = 1, q*w_2*y_4 = 1, a*y_3 = 1, -2*v_1*x_1 = 1, x_2*v_0 = 1, -20*w_0*x_3*y_3 = 1, v_1*x_1 = 1, -w_0*x_0*y_3 = 1, h*z_4 = 1, -w_2*x_0*y_0 = 1, q*w_2*y_0 = 1, -5*q*w_4*y_1 = 1, -12*w_2*x_1*y_1 = 1, -w_3 = 1, d*x_1 = 1, -x_4*v_0 = 1, -30*x_1*y_2*w_2 = 1, v_3 = 1, -w_0*x_0*y_0 = 1, -6*w_0*x_5*y_1 = 1, q*w_0*y_4 = 1, -3*w_0*x_2*y_1 = 1, -lm = 1, -6*w_2*x_0*y_2 = 1, q*w_5*y_0 = 1, -v_0*x_1 = 1, q*w_0*y_2 = 1, -6*v_2*x_2 = 1, w_2*y_1*q = 1, q*w_5*y_1 = 1, -w_0*x_0*y_4 = 1, w_2 = 1, -x_2*v_0 = 1, h*z_0 = 1, -k*y_4 = 1, -q*w_0*y_4 = 1, -10*q*w_3*y_2 = 1, q*w_0*y_1 = 1, w_3 = 1, -3*w_0*x_1*y_2 = 1, z_3 = 1, d*x_2 = 1, -30*w_2*x_2*y_1 = 1, -v_5*x_0 = 1, -60*w_2*x_1*y_3 = 1, y_2 = 1, q*w_0*y_5 = 1, -436608920299180809179630557104527219900824252088990974277221602192823366649534158 = 1, -q*w_0*y_1 = 1, w_6 = 1, -4*w_1*x_3*y_0 = 1, x_4 = 1, -w_0*x_2*y_0 = 1, v_1*x_3 = 1, -w_0*x_3*y_0 = 1, x_5 = 1, -w_5*x_0*y_0 = 1, -15*w_2*x_0*y_4 = 1, -w_0 = 1, d*x_0 = 1, -30*w_1*x_2*y_2 = 1, -6*w_0*x_1*y_5 = 1, q*w_1*y_5 = 1, v_2 = 1, -z_0 = 1, b*w_0 = 1, -10*w_0*x_2*y_3 = 1, x_2*v_1 = 1, x_3*v_0 = 1, y_3 = 1, -4*q*w_1*y_3 = 1, -20*w_1*x_3*y_1 = 1, -k*y_1 = 1, -12*w_1*x_2*y_1 = 1, -v_0*x_0 = 1, -z_4 = 1, -w_1*y_0*q = 1, -4*v_3*x_1 = 1, q*w_3*y_0 = 1, v_4*x_1 = 1, -z_3 = 1, q*w_1*y_3 = 1, -w_0*x_0*y_5 = 1, v_0*x_1 = 1, -5*w_4*x_1*y_0 = 1, -6*w_1*x_5*y_0 = 1, -2*y_1*w_1*x_0 = 1, q*w_0*y_0 = 1, -6*y_2*q*w_2 = 1, q*w_6*y_0 = 1, w_4 = 1, y_1*w_1*q = 1, q*w_1*y_4 = 1, -v_4*x_0 = 1, -15*w_0*x_2*y_4 = 1, -2*w_1*x_1*y_0 = 1, -k*y_3 = 1, -v_2*x_0 = 1, -12*w_1*x_1*y_2 = 1, -w_1 = 1, -3*w_1*x_0*y_2 = 1, -5*w_1*x_4*y_0 = 1, u*v_0 = 1, q*w_4*y_1 = 1, -30*w_1*x_4*y_1 = 1, -w_6 = 1, -q*w_2*y_0 = 1, -1 = 1, -5*w_0*x_1*y_4 = 1]
# 273, -5.548535779
# 1
# [c = 133, beta = 42]
# [c = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], beta = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]]