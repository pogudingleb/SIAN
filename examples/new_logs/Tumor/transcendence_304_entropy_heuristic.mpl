infolevel[Groebner]:=10;
Et_hat := [200353028956373442934091-x5_0, -k7*x1_0+x5_1, 63111623599406212019375-d_0, d_1, -x5_1+39173318450557368429281373248643490962474524227, -k7*x1_1+x5_2, k3*x1_0-k4*x2_0+k7*x1_0+x1_1, -x5_2-8292410455266503711298844313142733665648797749113070552169852521171927, -k7*x1_2+x5_3, (k3+k7)*x1_1+x1_2-k4*x2_1, b*d_0*k5*x2_0+a*k5*x2_0-k5*x2_0*x3_0-k5*x2_0*x4_0-k3*x1_0+k4*x2_0-k6*x3_0-k6*x4_0+x2_1, -x5_3-2519308262035224210550355375543242832043759994782638276681358979520466080683685761739887668710668050343196432499552295286030654514335794403, -k7*x1_3+x5_4, (k3+k7)*x1_2+x1_3-k4*x2_2, ((b*d_0+a-x3_0-x4_0)*x2_1-x2_0*(-b*d_1+x3_1+x4_1))*k5+k4*x2_1-k6*x3_1-k6*x4_1-x1_1*k3+x2_2, -a*k5*x2_0+k5*x2_0*x3_0+k6*x3_0+x3_1, -b*d_0*k5*x2_0+k5*x2_0*x4_0+k6*x4_0+x4_1, -x5_4+2456980888588886203469181703744780250555304935609793433199015879809525521994061271246247738891763885292136155234012751782040366159417590044019080686598889471157950231952497213348083295418557221377810596459618, -k7*x1_4+x5_5, (k3+k7)*x1_3+x1_4-k4*x2_3, ((b*d_0+a-x3_0-x4_0)*x2_2+(b*d_2-x3_2-x4_2)*x2_0-2*x2_1*(-b*d_1+x3_1+x4_1))*k5-x1_2*k3+k4*x2_2-k6*x3_2-k6*x4_2+x2_3, d_2, ((-a+x3_0)*x2_1+x2_0*x3_1)*k5+k6*x3_1+x3_2, ((-b*d_1+x4_1)*x2_0-x2_1*(b*d_0-x4_0))*k5+k6*x4_1+x4_2, -x5_5-2396195486619107806612855990228435485775996474815057645368178084175854650560873150384000436102337975179514985970417454555000818914522595638084689033466935874133911579943025105123713422371521218387913633459223372071229059606619568775001479608732041761025149048274787804545276905, -k7*x1_5+x5_6, (k3+k7)*x1_4+x1_5-k4*x2_4, ((b*d_0+a-x3_0-x4_0)*x2_3+(3*d_1*x2_2+3*d_2*x2_1+d_3*x2_0)*b+(-x3_3-x4_3)*x2_0+(-3*x3_2-3*x4_2)*x2_1-3*x2_2*(x3_1+x4_1))*k5-x1_3*k3+k4*x2_3-k6*x3_3-k6*x4_3+x2_4, d_3, ((-a+x3_0)*x2_2+x2_0*x3_2+2*x2_1*x3_1)*k5+k6*x3_2+x3_3, ((-d_0*x2_2-2*d_1*x2_1-d_2*x2_0)*b+x4_0*x2_2+2*x2_1*x4_1+x4_2*x2_0)*k5+k6*x4_2+x4_3, -x5_6+2336913907943107492502452345537947591083244358448586717024698762825381095242521826347694403462588688193224007267884461964141364149288115879702567850298957679868945743499845392302421976067096505589453481123161367675744240227736910914083919026741460426236944574408941055994267431648691758043300182484292809666336347730302545344508573590644559467100, -k7*x1_6+x5_7, (k3+k7)*x1_5+x1_6-k4*x2_5, ((b*d_0+a-x3_0-x4_0)*x2_4+(4*d_1*x2_3+6*d_2*x2_2+4*d_3*x2_1+d_4*x2_0)*b+(-x3_4-x4_4)*x2_0+(-4*x3_3-4*x4_3)*x2_1+(-6*x3_2-6*x4_2)*x2_2-4*x2_3*(x3_1+x4_1))*k5-x1_4*k3+k4*x2_4-k6*x3_4-k6*x4_4+x2_5, d_4, ((-a+x3_0)*x2_3+3*x3_2*x2_1+x3_3*x2_0+3*x2_2*x3_1)*k5+k6*x3_3+x3_4, ((-d_0*x2_3-3*d_1*x2_2-3*d_2*x2_1-d_3*x2_0)*b+3*x4_1*x2_2+x2_3*x4_0+3*x4_2*x2_1+x4_3*x2_0)*k5+k6*x4_3+x4_4, -x5_7-2279098948159406861136883335307685734204996206324585678672209217629104237337176202768305848266718635052537134833679716358140063022909869836314348026537162028957974369790760580153762786184369642329512476619229005265406504568235358934458832229678403602539174795241778501973355412571709550850196897322778780189343760193114034253211690107811935739684207599170602408542009049268686378383379926873049555352926680122510654, -k7*x1_7+x5_8, (k3+k7)*x1_6+x1_7-k4*x2_6, ((d_0*x2_5+5*d_1*x2_4+10*d_2*x2_3+10*d_3*x2_2+5*d_4*x2_1+d_5*x2_0)*b+(a-x3_0-x4_0)*x2_5+(-x3_5-x4_5)*x2_0+(-5*x3_4-5*x4_4)*x2_1+(-10*x3_3-10*x4_3)*x2_2+(-10*x3_2-10*x4_2)*x2_3-5*x2_4*(x3_1+x4_1))*k5-x1_5*k3+k4*x2_5-k6*x3_5-k6*x4_5+x2_6, d_5, ((-a+x3_0)*x2_4+4*x3_1*x2_3+6*x3_2*x2_2+4*x3_3*x2_1+x3_4*x2_0)*k5+k6*x3_4+x3_5, ((-d_0*x2_4-4*d_1*x2_3-6*d_2*x2_2-4*d_3*x2_1-d_4*x2_0)*b+4*x4_3*x2_1+6*x4_2*x2_2+4*x4_1*x2_3+x4_0*x2_4+x4_4*x2_0)*k5+k6*x4_4+x4_5, -x5_8+2222714323298798440107753255959597036998334023200490637352935676108665548254732314693142256670319737594437862443855490995027723474756207011566562193145020449103391825896872056340301235384524508473162800894268560318175464851039940769592407703703097005161691213418275814589829205988687004122292503892151055856872472986232391328143587058250895722192230484568779339450800625946148176002352101385132057158186742737159052942554219124989616543932080498140610146017487931008557686253786502142, -k7*x1_8+x5_9, (k3+k7)*x1_7+x1_8-k4*x2_7, ((d_0*x2_6+6*d_1*x2_5+15*d_2*x2_4+20*d_3*x2_3+15*d_4*x2_2+6*d_5*x2_1+d_6*x2_0)*b+(a-x3_0-x4_0)*x2_6+(-x3_6-x4_6)*x2_0+(-6*x3_5-6*x4_5)*x2_1+(-15*x3_4-15*x4_4)*x2_2+(-20*x3_3-20*x4_3)*x2_3+(-15*x3_2-15*x4_2)*x2_4-6*x2_5*(x3_1+x4_1))*k5-x1_6*k3+k4*x2_6-k6*x3_6-k6*x4_6+x2_7, d_6, ((-a+x3_0)*x2_5+5*x3_1*x2_4+10*x3_2*x2_3+10*x3_3*x2_2+5*x3_4*x2_1+x3_5*x2_0)*k5+k6*x3_5+x3_6, ((-d_0*x2_5-5*d_1*x2_4-10*d_2*x2_3-10*d_3*x2_2-5*d_4*x2_1-d_5*x2_0)*b+x2_0*x4_5+5*x4_4*x2_1+10*x4_3*x2_2+10*x4_2*x2_3+5*x4_1*x2_4+x4_0*x2_5)*k5+k6*x4_5+x4_6, -x5_9-2167724647052965640849258135430777385159800091058947245676507583810999206941665713929242799435840804540022596772571852443804587142424810883166330303768034349660463796396266320583769575259477618402341021181243943729485216660688941130937238063493192487951409539576579036065918840956640168359799062740262111946114471943526082013310527652473073226631585859498760174838318778008193428474697080378055036915550109917189610449770855150373826711632450227594679535134579277523153594844583807081528190312800651041530828730132645859393901142078389035270258980625992, -k7*x1_9+x5_10, (k3+k7)*x1_8+x1_9-k4*x2_8, ((d_0*x2_7+7*d_1*x2_6+21*d_2*x2_5+35*d_3*x2_4+35*d_4*x2_3+21*d_5*x2_2+7*d_6*x2_1+d_7*x2_0)*b+(a-x3_0-x4_0)*x2_7+(-x3_7-x4_7)*x2_0+(-7*x3_6-7*x4_6)*x2_1+(-21*x3_5-21*x4_5)*x2_2+(-35*x3_4-35*x4_4)*x2_3+(-35*x3_3-35*x4_3)*x2_4+(-21*x3_2-21*x4_2)*x2_5-7*x2_6*(x3_1+x4_1))*k5-x1_7*k3+k4*x2_7-k6*x3_7-k6*x4_7+x2_8, d_7, ((-a+x3_0)*x2_6+6*x3_1*x2_5+15*x3_2*x2_4+20*x3_3*x2_3+15*x3_4*x2_2+6*x2_1*x3_5+x3_6*x2_0)*k5+k6*x3_6+x3_7, ((-d_0*x2_6-6*d_1*x2_5-15*d_2*x2_4-20*d_3*x2_3-15*d_4*x2_2-6*d_5*x2_1-d_6*x2_0)*b+x2_0*x4_6+6*x4_5*x2_1+15*x4_4*x2_2+20*x4_3*x2_3+15*x4_2*x2_4+6*x2_5*x4_1+x4_0*x2_6)*k5+k6*x4_6+x4_7, -x5_10+2114095408566463828438376581226182659110362616909109475695611241008759000430321776114618831595630577636738132576506802879639332618048933204801836638404269022152744209246466685545042324872480036415233526115526124702998137674822488322617280042435349773753047452497305301526415584348408915141456852656550810836047037846495392170608137737756989329141353881018778912399068602551796018382230803649119027387200804569599055074052280119570913440091672292344072218751873485962302912720757391933420529392616830389153023823439024180735171194943529987919522111466056405058258696439066288891102714403343992704241637726058634236920198515, -d_1, -d_2, -d_3, -d_4, -d_5, -d_6, -d_7, 339214290049212301490536-x3_0, 201052406909511481727150839350642566397330463102402812262133724748744755090407886156158583790204179818775668141634433409831623736823446589829398293486995511064455198399740295313215607916021361202605672463513980124107929248249190100855401616023116322964496006471940562859221768468401107317847941989516274868540983170098294535738625167250604814718391292486683159326796978304031517308992611481305502510183344833740353734621341323351510616764636349468780407609690670481336831425877114027985068020778263512487644-x4_7, z_aux-1];
vars:=[x5_10, x5_9, x1_9, x5_8, x2_8, x1_8, x5_7, x4_7, x3_7, x2_7, x1_7, d_7, x5_6, x4_6, x3_6, x2_6, x1_6, d_6, x5_5, x4_5, x3_5, x2_5, x1_5, d_5, x5_4, x4_4, x3_4, x2_4, x1_4, d_4, x5_3, x4_3, x3_3, x2_3, x1_3, d_3, x5_2, x4_2, x3_2, x2_2, x1_2, d_2, x5_1, x4_1, x3_1, x2_1, x1_1, d_1, x5_0, x4_0, x3_0, x2_0, x1_0, d_0, z_aux, w_aux, a, b, k3, k4, k5, k6, k7];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [d, x3_0, x4_7] 115
# 297, [1 = 8, -k6 = 2, -k5*x2_0 = 2, k7*x1_2 = 1, -k4*x2_2 = 1, -k4*x2_7 = 1, -15*k5*x2_4*x3_2 = 1, k5*x2_3*x4_0 = 1, -k4*x2_8 = 1, -35*k5*x2_3*x4_4 = 1, -k6*x4_1 = 1, -k5*x2_0*x4_6 = 1, -a*k5*x2_5 = 1, k5*x2_1*x4_2 = 1, k5*x2_0*x4_3 = 1, -k4*x2_4 = 1, k4*x2_2 = 1, x5_5 = 1, -x1_4*k3 = 1, x5_2 = 1, -10*k5*x2_3*x4_2 = 1, k6*x3_5 = 1, x2_2*k5 = 1, x3_4 = 1, -k6*x3_4 = 1, x2_4*k5 = 1, -6*k5*x2_2*x3_2 = 1, -k7*x1_7 = 1, x1_5*k3 = 1, k5*x2_0*x4_4 = 1, -k5*x2_0*x3_5 = 1, -6*k5*x2_1*x3_5 = 1, x4_4 = 1, k5*x2_0*x4_2 = 1, -k6*x4_3 = 1, -35*k5*x2_3*x3_4 = 1, -a*k5*x2_0 = 1, -x2_7*k5 = 1, k6*x4_4 = 1, k7*x1_7 = 1, k6*x3_4 = 1, -x2_5*k5 = 1, -4*k5*x2_1*x3_3 = 1, k7*x1_4 = 1, -a*k5*x2_4 = 1, k5*x2_4*x4_1 = 1, x5_7 = 1, -k4*x2_0 = 1, k5*x2_4*x3_1 = 1, x1_4 = 1, -k6*x3_1 = 1, x1_2 = 1, k5*x2_1*x4_4 = 1, x1_1 = 1, -k5*x2_0*x3_6 = 1, -k4*x2_3 = 1, k4*x2_6 = 1, -15*k5*x2_4*x4_2 = 1, -k4*x2_6 = 1, x2_5*k5 = 1, -k5*x2_0*x3_3 = 1, -k6*x3_3 = 1, x3_3 = 1, -k7*x1_5 = 1, x5_3 = 1, -k7*x1_9 = 1, -x2_3*b*k5 = 1, -21*k5*x2_2*x4_5 = 1, k5*x2_3*x4_1 = 1, a*k5*x2_1 = 1, k5*x2_4*x3_2 = 1, -5*k5*x2_4*x3_1 = 1, x2_7 = 1, -6*k5*x2_5*x3_1 = 1, x1_3*k3 = 1, k5*x2_1*x3_2 = 1, -k7*x1_1 = 1, k5*x2_1*x4_5 = 1, -4*k5*x2_1*x4_3 = 1, -x2_0*b*k5 = 1, -2*k5*x2_1*x4_1 = 1, x2_8 = 1, -a*k5*x2_1 = 1, x1_7*k3 = 1, x2_6 = 1, -k6*x4_4 = 1, -k7*x1_6 = 1, -x2_2*b*k5 = 1, -x5_7 = 1, -x2_2*k5 = 1, x2_5*b*k5 = 1, -k5*x2_0*x4_4 = 1, k5*x2_2*x3_1 = 1, k5*x2_0*x3_3 = 1, x2_1 = 1, k5*x2_0*x3_4 = 1, k7*x1_6 = 1, x2_2*b*k5 = 1, x3_6 = 1, x1_8 = 1, k5*x2_1*x3_1 = 1, x1_5 = 1, -k7*x1_2 = 1, -6*k5*x2_5*x4_1 = 1, -k6*x4_2 = 1, -5*k5*x2_4*x4_1 = 1, x4_6 = 1, z_aux = 1, x5_9 = 1, x4_3 = 1, -k6*x3_2 = 1, x1_9 = 1, -k7*x1_8 = 1, -k5*x2_0*x4_3 = 1, x5_8 = 1, x2_3*k5 = 1, x2_1*k5 = 1, -k7*x1_3 = 1, x2_6*b*k5 = 1, x1_4*k3 = 1, -35*k5*x2_4*x4_3 = 1, -x1_1*k3 = 1, -35*k5*x2_4*x3_3 = 1, -k5*x2_2*x4_0 = 1, x2_6*k5 = 1, x3_1 = 1, x2_4 = 1, -6958898730416173093881179172160269835899249192578529656318533083801524084143529547565567783656836109381594164914365128709209750713714355577473815260343320536566315540021371260219115840733848346539882148721514624345490094875222521492456172344464943363725714853324139892775150896366671173178625873299254295613265111656970047959223811678076420145907653294534944695313803232019429145978205643400227003419669778944804877988674644 = 1, x3_2 = 1, k5*x2_3*x4_2 = 1, k7*x1_3 = 1, -x5_0 = 1, x2_7*b*k5 = 1, -k5*x2_5*x4_0 = 1, -x5_9 = 1, k5*x2_1*x3_5 = 1, -x5_5 = 1, a*k5*x2_5 = 1, -a*k5*x2_6 = 1, -3*k5*x2_1*x4_2 = 1, -k5*x2_1*x4_0 = 1, a*k5*x2_0 = 1, k5*x2_3*x4_3 = 1, -3*k5*x2_1*x3_2 = 1, k6*x4_0 = 1, k4*x2_4 = 1, -k4*x2_5 = 1, -x2_5*b*k5 = 1, -x5_2 = 1, k5*x2_1*x3_4 = 1, -x1_2*k3 = 1, -k5*x2_0*x3_4 = 1, a*k5*x2_6 = 1, x4_5 = 1, k4*x2_5 = 1, -x1_5*k3 = 1, k5*x2_2*x4_0 = 1, -x2_4*b*k5 = 1, -10*k5*x2_2*x4_3 = 1, k5*x2_2*x4_1 = 1, a*k5*x2_3 = 1, -k4*x2_1 = 1, -6*k5*x2_1*x4_5 = 1, -x2_1*b*k5 = 1, -x5_6 = 1, -k6*x3_6 = 1, k4*x2_3 = 1, -6*k5*x2_2*x4_2 = 1, k6*x3_6 = 1, k5*x2_4*x4_0 = 1, k5*x2_5*x3_1 = 1, -x1_7*k3 = 1, k7*x1_0 = 1, x1_7 = 1, k6*x4_3 = 1, k5*x2_2*x4_2 = 1, k5*x2_0*x3_1 = 1, k5*x2_3*x3_3 = 1, x3_7 = 1, -k5*x2_0*x3_7 = 1, x1_1*k3 = 1, -x2_1*k5 = 1, -3*k5*x2_2*x3_1 = 1, -7*k5*x2_6*x3_1 = 1, -k5*x2_0*x4_1 = 1, -10*k5*x2_3*x3_2 = 1, x4_1 = 1, x2_2 = 1, x2_0*b*k5 = 1, k5*x2_5*x4_0 = 1, -x1_6*k3 = 1, k6*x3_2 = 1, -k5*x2_0*x3_1 = 1, -3*k5*x2_2*x4_1 = 1, -x5_10 = 1, k5*x2_3*x3_2 = 1, -4*k5*x2_3*x4_1 = 1, k4*x2_7 = 1, -10*k5*x2_2*x3_3 = 1, -a*k5*x2_2 = 1, k6*x4_1 = 1, -21*k5*x2_5*x3_2 = 1, k5*x2_0*x4_5 = 1, -a*k5*x2_3 = 1, k5*x2_0*x3_6 = 1, -k7*x1_0 = 1, -k6*x3_5 = 1, k7*x1_1 = 1, k5*x2_0 = 1, -x2_6*b*k5 = 1, k5*x2_5*x4_1 = 1, -21*k5*x2_5*x4_2 = 1, -x2_4*k5 = 1, -x5_3 = 1, -x5_8 = 1, k5*x2_2*x3_2 = 1, -k6*x4_5 = 1, x1_6*k3 = 1, k5*x2_0*x4_1 = 1, -15*k5*x2_2*x3_4 = 1, x2_3 = 1, -k5*x2_6*x4_0 = 1, x5_1 = 1, a*k5*x2_4 = 1, a*k5*x2_2 = 1, k6*x4_6 = 1, -4*k5*x2_3*x3_1 = 1, -12663065265123807982888787825085092295220104869264970534300457573099384454871871971681419512550389673911239625647570704193591824310271170163587478897456228577628952185203849016366253155477248429754025782358995112199455497754634797335400427931744926077184451542373718207432386589563799194467605057247317745444113459220597257748557334597278589967366300139434366908915650091358478522218016413159746368505297248857162392719430084995529611140465692234972767094761025232312985037972811826503227149351368548411780454923678140192656211778714964364811316335115215560594285348 = 1, x1_2*k3 = 1, -x2_3*k5 = 1, x2_1*b*k5 = 1, -7*k5*x2_1*x3_6 = 1, x1_3 = 1, x5_10 = 1, k5*x2_0*x4_0 = 1, x2_5 = 1, k6*x4_2 = 1, -k6*x3_7 = 1, k5*x2_1*x3_3 = 1, k5*x2_1*x4_0 = 1, -20*k5*x2_3*x3_3 = 1, x5_4 = 1, a*k5*x2_7 = 1, k5*x2_0*x3_2 = 1, k6*x3_1 = 1, x2_3*b*k5 = 1, -2*k5*x2_1*x3_1 = 1, k5*x2_6*x4_0 = 1, k6*x4_5 = 1, x2_4*b*k5 = 1, -x5_4 = 1, -15*k5*x2_2*x4_4 = 1, -x5_1 = 1, k5*x2_4*x4_2 = 1, k7*x1_8 = 1, -2101569993595117272019692832471287042789533578938642789163863463052654610834974025134185456315817485350439409513875194730849959472086459311188 = 1, k5*x2_3*x3_1 = 1, k7*x1_5 = 1, -k5*x2_0*x4_2 = 1, k5*x2_2*x4_4 = 1, k5*x2_2*x3_3 = 1, -7*k5*x2_1*x4_6 = 1, -k3*x1_0 = 1, -k5*x2_3*x4_0 = 1, k5*x2_0*x4_6 = 1, x1_8*k3 = 1, k5*x2_0*x3_5 = 1, -k5*x2_7*x4_0 = 1, k3*x1_0 = 1, -5*k5*x2_1*x4_4 = 1, k4*x2_0 = 1, -20*k5*x2_3*x4_3 = 1, k4*x2_1 = 1, k5*x2_2*x4_3 = 1, -k7*x1_4 = 1, k5*x2_2*x3_4 = 1, x3_5 = 1, -k5*x2_0*x4_5 = 1, -7*k5*x2_6*x4_1 = 1, -3824214005558499939308261991650147977828330559559408503020592143973822885742745186147450083465592143287134790955972080835225780889448170397953102928712428502720686435894501620622866189839031176439899083141240277971985493390052575656538797740462078805039725342477858806823960469070084 = 1, -k6*x4_0 = 1, -x2_6*k5 = 1, -5*k5*x2_1*x3_4 = 1, -21*k5*x2_2*x3_5 = 1, k6*x3_3 = 1, k5*x2_1*x4_3 = 1, -x1_3*k3 = 1, x4_2 = 1, -k5*x2_0*x3_2 = 1, x1_6 = 1, k6 = 1, -k6*x4_6 = 1, -k5*x2_0*x4_0 = 1, -k5*x2_4*x4_0 = 1, x5_6 = 1, -1 = 1, k5*x2_1*x4_1 = 1]
# 306, -5.660159869
# 3
# [d = 15, x3_0 = 17, x4_7 = 3]
# [d = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], x3_0 = [3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3], x4_7 = [3, 2, 1]]