infolevel[Groebner]:=10;
Et_hat := [6752603289226891498-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-405365731559319644452568698063512552512587609299847603839/9556902864018205328, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+9267357444488031986958665921521027381147425517484071965246135657519107360857614141133911619803187236622325572635/34783031502679892228472151357002023264528910116643988736, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-105933860853923500618613276159915121892189096115215395575546981326408990319530991264870186567038743898347748991848252444944720544031241390367116295520912401372878770703/63297665453499009848461689929290398002461660650803159137859613888934173105532134149686921216, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0)*gama)*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+12026405480481703345387030871551648378599049711117726816688880801179022456157026549006311231629909314787483403320172197535158599519326538595303574008574471134568956479560994965058743080501183612964513449561259846629558118901/1151881902977430266855813519007328245886272164374384990963916555034094702359719209569985034588685730292913017379166784110505000960, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-69310313650001498750851733760335659510685714421550624344643881880387512916728570558375307067952222196045710794137817678789933620641705751441687877082451849044480903224270239429678637706451035178964588528697844924900266662098239867868161085558041473023647596851670735990727009043675639692976679547/5240445713521811214066218526633723168023560730883718974490763250526750021322094197258839896422044692454730006020097603405755466711038842286590842111817165016376934400, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+42153872377235121330870058087409055158319130747070994031998134015673846209866664174509116488498499113851755799423621515403322589775671341721190417892923000321591262031665013700281409077552701847372423026691711724769526657212232817746827089655835856114863644788811027115091439412413847273246564795310720004417612455539605913863447299739372375234555336184996178348738316673/95364884908369690915081795268995478407075227281945896893008248321807361789167010606003877756502484811072750051552674868410212741844046169077203846320907852102912074722271994738419169501395763326091264000, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+535356589308077377571310676591788216338028835422705918268863027726009154455647163162967152308611125666305205717305212974333463288194256244796237845170492884058784064378670925844347413023566425232125382748815247698906805157727737533900778582171377442831882629194179902430102314906465954319955802044313639848572476017674504947632985308698056075211321570581542378610169800859594084401884364425292512904140313509651700739352112034275173284866691521/173543659657046941614017077573018380944123000752245126011636267992868252194359965774234294953545371398735176545198811000433069684408601487996633007386446240995700351362316629312860029395326294048169440327359539253659261518180691197558784000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-3050802491428429674557835478588920340696767149043348800125695069286734178803180765873505572111273443419634944379629292680255925359823813465844981227128032543617142453377454973980869901059339188757091296414458356329293449015411188243553495840858497900026724491592239121153231665889899258661715219081960413435330913493546521968103857627433713509145544395814192390032952867465699123421236157887841088866559742938845968905363768497112210816198708682327141589705917264916027964796033975177667041161821027426506081752708159751/3158122807582572882598720030344479700482218996391813882887124921303715167422199891891889018753769063947955968183019197531072843935462396259763174691722793187462893399775309534744057440176282166426748924990114112553360255898262594826978444412809996460199755959861018403799040000, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+127368217757609347810358306485793659511127141816851795175982920839264047162590798330308940377601158739361086937319787293099759832594401386440940107771463375851740275726141846562504382965630765051375658341045258632847576797990718944616642149930191796370472077990297755210158403154838127803385163985599440390759662992131088278031482267489932955123904874737227509187312771282659213359399319263295844804317725499933850780226009565484752863111729062823113356275171890653553123608105439422391412295808271758554149897997785772560638044339322218427224596666331739789877343039459885593236711442236335867/3591941822983878957202470072507957268498192701616906664474641041638132650160201648316956450761750476551224850849252527429819713821532451831607127905086722680760817369451737914349525925959499571735077162749476441049856204593546290617098927661526828783056527013591248433559900926331574669662196948658433713766400000, 6663190258058088688795669409377249239-x2_1, 9009904715707953094731448979707470147-x3_1, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_1, x3_1] 78
# 274, [1 = 8, -1 = 2, delta*x3_2 = 1, x2_2 = 1, x2_5 = 1, x1_3*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, x4_5*delta*sgm = 1, b*x1_0*x4_5 = 1, -gama = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, delta*x3_0 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, delta*sgm*x3_0*x4_0 = 1, -gama*sgm*x2_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, x4_6*delta*sgm = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x3_0*x4_3 = 1, -gama*x2_2 = 1, -21*gama*sgm*x2_5*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x3_4 = 1, delta*sgm*x3_4*x4_1 = 1, -2*gama*x4_1*sgm = 1, b*c*x1_7 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_1*x4_3 = 1, -7*gama*x4_6*sgm = 1, x4_6 = 1, x3_5*x4_2 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, x3_0*x4_5 = 1, b*x1_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_7*x4_0 = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -gama*x2_6 = 1, x1_3*x4_5 = 1, -gama*sgm*x2_0*x4_3 = 1, -alpha*x1_2 = 1, x1_6*x4_0 = 1, -x1_4 = 1, -gama*x2_0 = 1, x1_1*c = 1, x1_1*x4_6 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, -4*gama*x4_3*sgm = 1, x4_3*delta*sgm = 1, delta*x4_0*sgm = 1, -gama*x2_5 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_0*x4_6 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, x3_7 = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, x4_7 = 1, -gama*sgm*x2_7*x4_0 = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, delta*sgm*x3_4*x4_3 = 1, x4_1 = 1, delta*sgm*x3_0*x4_7 = 1, x1_1*x4_7 = 1, x3_0*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, b*x1_3*x4_0 = 1, beta*x2_5 = 1, delta*sgm*x3_3*x4_2 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, -4*gama*sgm*x2_3*x4_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_3*x4_2 = 1, x1_3*x4_4 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_4*x4_0 = 1, -x1_0 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, x2_3 = 1, x1_9*x4_0 = 1, -15*gama*sgm*x2_2*x4_4 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_4*x4_0 = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, -20*gama*sgm*x2_3*x4_3 = 1, x3_2*x4_3 = 1, beta*x2_4 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, delta = 1, x2_4 = 1, delta*x3_5 = 1, x1_5*x4_4 = 1, -gama*sgm*x2_0*x4_7 = 1, x1_5*x4_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_0*x4_2 = 1, -alpha*x1_5 = 1, delta*sgm*x3_5*x4_1 = 1, x3_3*x4_5 = 1, b*c*x1_5 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -gama*sgm*x4_0 = 1, -6*gama*sgm*x2_2*x4_2 = 1, b*x1_0*x4_7 = 1, b*x1_1*x4_2 = 1, x3_0*x4_6 = 1, x1_3*c = 1, x1_4*x4_2 = 1, x1_5*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x3_3 = 1, beta*x2_0 = 1, x3_0*x4_8 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, x3_0*x4_1 = 1, z_aux*x3_0*c = 1, beta = 1, -5*gama*x4_4*sgm = 1, b*x1_5*x4_1 = 1, x4_2 = 1, b*x1_1*x4_7 = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, x4_4 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, x1_3*x4_0 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_1*x4_1 = 1, x3_5*x4_1 = 1, b*x1_1*x4_4 = 1, x3_3*x4_1 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_0*x4_6 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_5*c = 1, x1_6*x4_2 = 1, -10*gama*sgm*x2_3*x4_2 = 1, b*x1_1*x4_5 = 1, b*c*x1_1 = 1, x3_5 = 1, b*c*x1_6 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, x4_3 = 1, -6*gama*x4_5*sgm = 1, x1_8*x4_1 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, b*x1_4*x4_1 = 1, delta*x4_1*sgm = 1, -3*x4_2*gama*sgm = 1, b*x1_2*x4_1 = 1, -alpha*x1_3 = 1, x3_7*x4_1 = 1, -gama*x2_4 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_4*x4_0 = 1, -gama*sgm*x2_2*x4_0 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_5 = 1, x1_4*x4_0 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, x4_2*delta*sgm = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, x4_4*delta*sgm = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, delta*sgm*x3_5*x4_0 = 1, x3_0*x4_4 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x4_5 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, x3_5*x4_3 = 1, x1_7*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, x3_0*x4_2 = 1, -gama*x2_3 = 1, x3_6*x4_2 = 1, x3_2 = 1]
# 282, -5.577999884
# 2
# [x2_1 = 10, x3_1 = 16]
# [x2_1 = [4, 1, 4, 2, 2, 4, 4, 4, 4, 4], x3_1 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2]]