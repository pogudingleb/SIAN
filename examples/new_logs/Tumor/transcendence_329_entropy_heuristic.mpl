infolevel[Groebner]:=10;
Et_hat := [137945467401764746507891-x5_0, -k7*x1_0+x5_1, 131187296234536167976960-d_0, d_1, -x5_1+195803424280433871444250597724363895958866561084, -k7*x1_1+x5_2, k3*x1_0-k4*x2_0+k7*x1_0+x1_1, -x5_2-60108634234935180697278854934156623875397986324755799006172298749086200, -k7*x1_2+x5_3, (k3+k7)*x1_1+x1_2-k4*x2_1, b*d_0*k5*x2_0+a*k5*x2_0-k5*x2_0*x3_0-k5*x2_0*x4_0-k3*x1_0+k4*x2_0-k6*x3_0-k6*x4_0+x2_1, -x5_3-226209605328005836002050442273276968328098070692984997788233940905952112785826556225592022998969527908548864810563023701384621593455732462816, -k7*x1_3+x5_4, (k3+k7)*x1_2+x1_3-k4*x2_2, ((b*d_0+a-x3_0-x4_0)*x2_1-x2_0*(-b*d_1+x3_1+x4_1))*k5+k4*x2_1-k6*x3_1-k6*x4_1-x1_1*k3+x2_2, -a*k5*x2_0+k5*x2_0*x3_0+k6*x3_0+x3_1, -b*d_0*k5*x2_0+k5*x2_0*x4_0+k6*x4_0+x4_1, -x5_4+5115346883522788421941175122816644969355945832973372230204118210546390679692976593567569548323463374686058501103825571293088152686241697893147123186018945269130878096665820482764157245386141572093288095993775280, -k7*x1_4+x5_5, (k3+k7)*x1_3+x1_4-k4*x2_3, ((b*d_0+a-x3_0-x4_0)*x2_2+(b*d_2-x3_2-x4_2)*x2_0-2*x2_1*(-b*d_1+x3_1+x4_1))*k5-x1_2*k3+k4*x2_2-k6*x3_2-k6*x4_2+x2_3, d_2, ((-a+x3_0)*x2_1+x2_0*x3_1)*k5+k6*x3_1+x3_2, ((-b*d_1+x4_1)*x2_0-x2_1*(b*d_0-x4_0))*k5+k6*x4_1+x4_2, -x5_5-115674901164449941361526762753429385250784705413731493489024090090991962518897721828577922101170191041427412793496913963413742199260576061053534219968689954697519875299738516595457745177407959182129832563135471174565138813203653196176985481669299258019241055614715225663568968977584, -k7*x1_5+x5_6, (k3+k7)*x1_4+x1_5-k4*x2_4, ((b*d_0+a-x3_0-x4_0)*x2_3+(3*d_1*x2_2+3*d_2*x2_1+d_3*x2_0)*b+(-x3_3-x4_3)*x2_0+(-3*x3_2-3*x4_2)*x2_1-3*x2_2*(x3_1+x4_1))*k5-x1_3*k3+k4*x2_3-k6*x3_3-k6*x4_3+x2_4, d_3, ((-a+x3_0)*x2_2+x2_0*x3_2+2*x2_1*x3_1)*k5+k6*x3_2+x3_3, ((-d_0*x2_2-2*d_1*x2_1-d_2*x2_0)*b+x4_0*x2_2+2*x2_1*x4_1+x4_2*x2_0)*k5+k6*x4_2+x4_3, -x5_6+2615791863989951177518799462190837413972783925256323991080701930367782218603680310906389831736363151222829703508359297961740473424645060931456148136246361471127004425578070264218796866361621831002368107525558718222259229913310980630947015112857019108239355291976728375171580867768383310972718276767339679077245471192644868098630925691080848233016112752, -k7*x1_6+x5_7, (k3+k7)*x1_5+x1_6-k4*x2_5, ((b*d_0+a-x3_0-x4_0)*x2_4+(4*d_1*x2_3+6*d_2*x2_2+4*d_3*x2_1+d_4*x2_0)*b+(-x3_4-x4_4)*x2_0+(-4*x3_3-4*x4_3)*x2_1+(-6*x3_2-6*x4_2)*x2_2-4*x2_3*(x3_1+x4_1))*k5-x1_4*k3+k4*x2_4-k6*x3_4-k6*x4_4+x2_5, d_4, ((-a+x3_0)*x2_3+3*x3_2*x2_1+x3_3*x2_0+3*x2_2*x3_1)*k5+k6*x3_3+x3_4, ((-d_0*x2_3-3*d_1*x2_2-3*d_2*x2_1-d_3*x2_0)*b+3*x4_1*x2_2+x2_3*x4_0+3*x4_2*x2_1+x4_3*x2_0)*k5+k6*x4_3+x4_4, -x5_7-59151700211859525244350282051754544908996423510505958126177149824591291316373871089685534430053263892798853430835884848208808190455640939654196028763274755739486020808424006888778261029484602938480222777197546804402576153353142878334374090886714733635283795245816708491477210737118008365703582180254870429313470878185983848286286881904447825765535216624884416772707052726443221869195793459854648167984244739587005024882448, -k7*x1_7+x5_8, (k3+k7)*x1_6+x1_7-k4*x2_6, ((d_0*x2_5+5*d_1*x2_4+10*d_2*x2_3+10*d_3*x2_2+5*d_4*x2_1+d_5*x2_0)*b+(a-x3_0-x4_0)*x2_5+(-x3_5-x4_5)*x2_0+(-5*x3_4-5*x4_4)*x2_1+(-10*x3_3-10*x4_3)*x2_2+(-10*x3_2-10*x4_2)*x2_3-5*x2_4*(x3_1+x4_1))*k5-x1_5*k3+k4*x2_5-k6*x3_5-k6*x4_5+x2_6, d_5, ((-a+x3_0)*x2_4+4*x3_1*x2_3+6*x3_2*x2_2+4*x3_3*x2_1+x3_4*x2_0)*k5+k6*x3_4+x3_5, ((-d_0*x2_4-4*d_1*x2_3-6*d_2*x2_2-4*d_3*x2_1-d_4*x2_0)*b+4*x4_3*x2_1+6*x4_2*x2_2+4*x4_1*x2_3+x4_0*x2_4+x4_4*x2_0)*k5+k6*x4_4+x4_5, -x5_8+1337615460205875028283553192294162035829166791502158463598707464715595102709246582797655705115585105626194511845260577805801960234708185697888996673966266653048782836404212590616052820370848411203242520321789420241897426817672930554190767818596676807454972822933581099460679010105229129287548045346204418592277472283817080267687699947977494450215314319684010389260787733958950414201011541134325950687272046372536275501253510010045993041678446275012716760246943730605537204317540661913728319760, -k7*x1_8+x5_9, (k3+k7)*x1_7+x1_8-k4*x2_7, ((d_0*x2_6+6*d_1*x2_5+15*d_2*x2_4+20*d_3*x2_3+15*d_4*x2_2+6*d_5*x2_1+d_6*x2_0)*b+(a-x3_0-x4_0)*x2_6+(-x3_6-x4_6)*x2_0+(-6*x3_5-6*x4_5)*x2_1+(-15*x3_4-15*x4_4)*x2_2+(-20*x3_3-20*x4_3)*x2_3+(-15*x3_2-15*x4_2)*x2_4-6*x2_5*(x3_1+x4_1))*k5-x1_6*k3+k4*x2_6-k6*x3_6-k6*x4_6+x2_7, d_6, ((-a+x3_0)*x2_5+5*x3_1*x2_4+10*x3_2*x2_3+10*x3_3*x2_2+5*x3_4*x2_1+x3_5*x2_0)*k5+k6*x3_5+x3_6, ((-d_0*x2_5-5*d_1*x2_4-10*d_2*x2_3-10*d_3*x2_2-5*d_4*x2_1-d_5*x2_0)*b+x2_0*x4_5+5*x4_4*x2_1+10*x4_3*x2_2+10*x4_2*x2_3+5*x4_1*x2_4+x4_0*x2_5)*k5+k6*x4_5+x4_6, -x5_9-30247906872895751967156559055799426531608907882060030760351075562381203228191505643358032431385249591269685099724914174041613648965285579131710994600410887944348861763704740427848362574083684675297406039486923131245093614184569582370194228445037971740941323279481362130423339095283811919439730257353989790598351035568409178275514109599316798813989188881828986227673915560404043421481079574571803840801399139041083492399338701509468768193931402054674520653178063012914901394342363208474665520572866680959456620827540677898529627307080556054985812824528984566584560, -k7*x1_9+x5_10, (k3+k7)*x1_8+x1_9-k4*x2_8, ((d_0*x2_7+7*d_1*x2_6+21*d_2*x2_5+35*d_3*x2_4+35*d_4*x2_3+21*d_5*x2_2+7*d_6*x2_1+d_7*x2_0)*b+(a-x3_0-x4_0)*x2_7+(-x3_7-x4_7)*x2_0+(-7*x3_6-7*x4_6)*x2_1+(-21*x3_5-21*x4_5)*x2_2+(-35*x3_4-35*x4_4)*x2_3+(-35*x3_3-35*x4_3)*x2_4+(-21*x3_2-21*x4_2)*x2_5-7*x2_6*(x3_1+x4_1))*k5-x1_7*k3+k4*x2_7-k6*x3_7-k6*x4_7+x2_8, d_7, ((-a+x3_0)*x2_6+6*x3_1*x2_5+15*x3_2*x2_4+20*x3_3*x2_3+15*x3_4*x2_2+6*x2_1*x3_5+x3_6*x2_0)*k5+k6*x3_6+x3_7, ((-d_0*x2_6-6*d_1*x2_5-15*d_2*x2_4-20*d_3*x2_3-15*d_4*x2_2-6*d_5*x2_1-d_6*x2_0)*b+x2_0*x4_6+6*x4_5*x2_1+15*x4_4*x2_2+20*x4_3*x2_3+15*x4_2*x2_4+6*x2_5*x4_1+x4_0*x2_6)*k5+k6*x4_6+x4_7, -x5_10+684005155002136778472057849060637947464186900094557119847429815576871583104624589057967290730891548484621714409910050869054018455562017387014196754677152921880464747288953494102740047434004124884334521641117753465466063570926636458184400968686425912282473224543423688365565385740160762435819409904145429984174325436651623933264508183869071605824927357160146902332650864607865146276765033144407251005683489667510841401903416736876650056987485420880009158162574138826211782066308599050468592197059092645661061065309228931131779720624986390512273795458587548530164547661149464304736848138827792979860058294392729727443081091991233858672, -d_1, -d_2, -d_3, -d_4, -d_5, -d_6, -d_7, -77813450296173872910313171426156514651716395324335229288378869468427234209681815370718411143649854535903370832334649596355505298965207859355-x3_2, 1562398255595314554426529352976578360631444370576026808115643956986331533571737881244133768540311720485540142273464811114462212321702817989955354172758422045629205670370465435535785831172538600006279285392277097122106579825383165128048758232991300973618970005787848086231482133777070578066248116795043621306679501622026749117394353693643199507619920031104127584619561757807574-x4_5, z_aux-1];
vars:=[x5_10, x5_9, x1_9, x5_8, x2_8, x1_8, x5_7, x4_7, x3_7, x2_7, x1_7, d_7, x5_6, x4_6, x3_6, x2_6, x1_6, d_6, x5_5, x4_5, x3_5, x2_5, x1_5, d_5, x5_4, x4_4, x3_4, x2_4, x1_4, d_4, x5_3, x4_3, x3_3, x2_3, x1_3, d_3, x5_2, x4_2, x3_2, x2_2, x1_2, d_2, x5_1, x4_1, x3_1, x2_1, x1_1, d_1, x5_0, x4_0, x3_0, x2_0, x1_0, d_0, z_aux, w_aux, a, b, k3, k4, k5, k6, k7];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [d, x3_2, x4_5] 118
# 170, [1 = 14, k5*x2_0 = 4, x2_1*k5 = 4, k6 = 4, x2_2*k5 = 3, k7*x1_5 = 2, k5*x2_2*x3_3 = 2, k6*x4_6 = 2, x5_10 = 2, b*x2_1*k5 = 2, a*k5*x2_1 = 2, a*k5*x2_4 = 2, k6*x4_3 = 2, k6*x3_4 = 2, x5_9 = 2, k6*x3_6 = 2, k4*x2_5 = 2, k5*x2_2*x3_1 = 2, x2_3*k5 = 2, k5*x2_0*x3_3 = 2, x2_5*b*k5 = 2, k5*x2_0*x4_1 = 2, k6*x4_0 = 2, k5*x2_1*x3_5 = 2, x5_1 = 2, k5*x2_4*x4_0 = 2, x2_4*k5 = 2, k6*x3_3 = 2, k5*x2_2*x4_2 = 2, x1_5*k3 = 2, k4*x2_6 = 2, x1_6*k3 = 2, b*k5*x2_0 = 2, x5_4 = 2, k4*x2_7 = 2, k6*x4_1 = 2, k7*x1_6 = 2, k5*x2_5*x3_0 = 2, k5*x2_0*x3_5 = 2, k6*x4_4 = 2, k5*x2_0*x3_0 = 2, k5*x2_2*x3_4 = 2, b*x2_4*k5 = 2, x1_1*k3 = 2, k5*x2_1*x4_3 = 2, b*x2_3*k5 = 2, k5*x2_6*x3_0 = 2, k7*x1_2 = 2, k4*x2_2 = 2, k5*x2_3*x3_0 = 2, x5_6 = 2, k5*x2_2*x4_1 = 2, k5*x2_6*x4_0 = 2, k5*x2_0*x4_3 = 2, k4*x2_3 = 2, k5*x2_4*x4_2 = 2, k5*x2_1*x3_3 = 2, k6*x3_1 = 2, k5*x2_1*x3_0 = 2, k5*x2_2*x4_4 = 2, k5*x2_0*x4_2 = 2, k5*x2_0*x4_6 = 2, a*k5*x2_5 = 2, x2_6*b*k5 = 2, k5*x2_5*x4_0 = 2, k5*x2_3*x3_1 = 2, k5*x2_3*x4_2 = 2, k5*x2_1*x4_4 = 2, k6*x3_0 = 2, x5_8 = 2, k5*x2_3*x4_1 = 2, k5*x2_0*x4_0 = 2, k7*x1_7 = 2, k5*x2_0*x3_6 = 2, b*x2_2*k5 = 2, k4*x2_0 = 2, k5*x2_2*x3_0 = 2, k5*x2_4*x4_1 = 2, k7*x1_3 = 2, k5*x2_0*x3_1 = 2, k5*x2_2*x4_0 = 2, k5*x2_1*x4_0 = 2, k5*x2_1*x3_1 = 2, x1_4*k3 = 2, k4*x2_1 = 2, x1_1*k7 = 2, k3*x1_0 = 2, k6*x3_5 = 2, x5_2 = 2, x5_7 = 2, k5*x2_5*x3_1 = 2, k7*x1_0 = 2, k5*x2_0*x4_4 = 2, k5*x2_3*x3_3 = 2, k6*x4_2 = 2, k4*x2_4 = 2, k5*x2_4*x3_1 = 2, x5_5 = 2, x5_3 = 2, k5*x2_1*x3_4 = 2, k5*x2_4*x3_0 = 2, k7*x1_8 = 2, a*k5*x2_6 = 2, k5*x2_2*x4_3 = 2, k5*x2_1*x4_2 = 2, x1_2*k3 = 2, a*k5*x2_3 = 2, k5*x2_0*x3_4 = 2, a*k5*x2_2 = 2, k5*x2_5*x4_1 = 2, k5*x2_3*x4_0 = 2, k5*x2_1*x4_1 = 2, x1_7*k3 = 2, k5*x2_3*x4_3 = 2, k7*x1_4 = 2, a*k5*x2_0 = 2, x1_3*k3 = 2, x3_4 = 1, x2_4 = 1, x2_1*k5*x4_6 = 1, x4_1 = 1, x2_6*k5*x4_1 = 1, a*k5*x2_7 = 1, x2_6*k5*x3_1 = 1, k5*x2_0*x4_7 = 1, x1_3 = 1, x3_7 = 1, z_aux = 1, x5_0 = 1, x2_2 = 1, x4_6 = 1, x2_8 = 1, x1_9 = 1, x4_2 = 1, x2_2*k5*x3_5 = 1, x3_6 = 1, k6*x3_7 = 1, x4_4 = 1, x2_5*k5*x4_2 = 1, x1_2 = 1, x3_1 = 1, x2_6 = 1, k5*x2_7*x3_0 = 1, k7*x1_9 = 1, x4_3 = 1, x2_3*k5*x3_4 = 1, k5*x2_7*x4_0 = 1, k4*x2_8 = 1, x4_7 = 1, x2_1 = 1, x2_4*k5*x3_3 = 1, x2_4*k5*x4_3 = 1, x3_3 = 1, x2_3 = 1, x3_5 = 1, x2_3*k5*x4_4 = 1, x2_5 = 1, x1_1 = 1, x2_5*k5 = 1, x1_8*k3 = 1, k6*x4_7 = 1, x2_1*k5*x3_6 = 1, x2_7 = 1, x1_4 = 1, k5*x2_0*x3_7 = 1, x1_7 = 1, x1_6 = 1, x1_5 = 1, x2_7*b*k5 = 1, x1_8 = 1]
# 306, -5.030306866
# 117
# [d = 15, x3_2 = 14, x4_5 = 8]
# [d = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], x3_2 = [3, 2, 1, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3], x4_5 = [3, 2, 1, 3, 3, 2, 3, 3]]