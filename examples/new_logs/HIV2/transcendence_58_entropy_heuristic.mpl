infolevel[Groebner]:=10;
Et_hat := [197838202819281962147-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 331432493648306562263-z_0, -c*q*w_0*y_0+h*z_0+z_1, 264546817835894454204-k_0, k_1, -w_1+821744351458497898034294174872244717952085495336716100449120454971902463471492928, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 186649071782291994186068358037050923089460107346781630926700768967215313094216654-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+3413212258962827006351710369423653995202923182858039260824731514091491159568247002374566747000635434165682986472306272924060142970298952645700, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, 775268973617675059745598830999656061238012815318065578464404068461651106332102360686814207876089905381417339651413337090961606796825904981196-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+14177180413904563642002169229126080129515505241422951350604325313664506951559809882057612886406184473704192026378471184251944813356556065426389265526323505489210594856472128577447869148878577693354880180, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 3220171285690937787263561203746791900602663048418957519514364358595971994895449943649833201532255535289658167180697410337183283340082151959695774672998148157594005880813194959546011606900128706051231832-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+58886593987997324091906197508019733241135905402590235548195633036958083573404049941438197952552012117385917275599539001158891356503645819911025649514732908640996545674750266153475293670267698626359185041253430835132798177646191224316881830854267647720161419814708, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, 13375361922199354927116539556026295317596844591352661029292848675079928860331766065183312419248965922240468246801255713004084613678079525582663666426315173936863760808274794154270652507164291934308888265413340988250135690601985386123754316006317754520654414230992-z_4, -4*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*q*c+h*z_4+z_5, -w_5+244592426016268484530619137556817744813774472917169763343364136731447836255549871697778100167036183357844263357373671200256206947938914166857610285130873224886537810285657703807642105792317427636147151525860420639708329139152856984000491296940649750739158426205025347160449440881814315449610802577360733868852201678219298996, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 55556146141907660046180284504701897318233696871919958096454668493691148646733933558927378309863913203650912785157866401310215539388626319215031409305262999971640686327067493455947637199946176838931169064137488367076019315937819639906455358557322023215228163640113832297660616116920624074476805737677013407059958542524827264-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+1015943541864856596590002228401339739622693112171500038048043015418571558820750285453402372146557010068727156705487598064373030589572933597284832791945588568501773288604060272047258466298813967206971357371358036855560711840596796574585623627076371115427880393307765360604829022495045690950772385696738114100426373207372800978181450700451899025069526167152956211453852362283334611124532, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7+4219841542387985477884714636725351990107480734860149528479477601366304931129452732584215595166865155088918814275678364618996876953296823252473807550533940889018585718084912623950480012465100951338131262272393306933066791522373393621163168828305783413773753987522950099963229441834760332457267974765818656606178567018194777973590216107022739923181682424025871319232598693845324786282157729306852178034078129702995241217517778431641670267060103220, 230759017370460835978494339505126733366319332659005616922662777601845336866286748485528954902439344273663163638790143342413185772376815954771560211736282102150299628281132726290980947900458312367022787898983102709784022140605409422266246371547225834135591154902673874528026345037647267631410371773083840967198320070423460480977553789014152797187840985282607016200341356445759451561568-z_6, -k_1, -k_2, -k_3, -k_4, 3628164716476494960322583115383313409833609693564718165403890734871390911041988734556621511792028351-x_2, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [k, x_2] 97
# 265, [1 = 9, -12*c*w_1*x_1*y_2 = 1, b*w_2 = 1, -y_0 = 1, a*y_0 = 1, -z_6 = 1, -2*c*w_1*x_1*y_0 = 1, -10*v_3*beta = 1, -5*c*w_4*x_0*y_1 = 1, u*v_2 = 1, -beta*v_3*x_0 = 1, beta*v_3*x_0 = 1, -60*c*w_2*x_3*y_1 = 1, beta*v_1*x_3 = 1, c*q*w_3*y_2 = 1, -3*c*w_0*x_1*y_2 = 1, -6*c*w_1*x_0*y_5 = 1, c*q*w_0*y_1 = 1, -beta*v_0*x_5 = 1, -6*v_2*beta = 1, w_1 = 1, h*z_2 = 1, beta*v_2*x_3 = 1, u*v_3 = 1, -30*c*w_1*y_2 = 1, x_6 = 1, v_2*beta = 1, d*x_3 = 1, -164721582855644909701359059267008809550018413877809439532012416060209601998657816842945745567897204576735231792678462717632891988364645853664005984739977361895622179151588377402130131036243839382931271497947324584784195495108097455350355250611728672076597977632656 = 1, -3*beta*v_1 = 1, -15*c*w_4*y_0 = 1, -4*beta*v_1*x_3 = 1, -c*q*w_0*y_1 = 1, z_6 = 1, -15*c*w_2*x_4*y_0 = 1, c*q*w_0*y_3 = 1, v_5 = 1, -6*c*w_0*x_5*y_1 = 1, -3*beta*v_2*x_1 = 1, c*q*w_3*y_3 = 1, c*q*w_0*y_4 = 1, -4*c*w_1*x_0*y_3 = 1, -c*w_0*x_0*y_5 = 1, b*w_6 = 1, z_1 = 1, c*q*w_0*y_6 = 1, -6*c*w_2*y_0 = 1, -c*w_0*x_4*y_0 = 1, u*v_4 = 1, v_4 = 1, -10*c*w_0*y_3 = 1, -20*c*w_3*x_0*y_3 = 1, beta*v_1*x_4 = 1, a*y_2 = 1, -beta*v_0*x_4 = 1, c*q*w_3*y_0 = 1, -c*w_2*x_0*y_0 = 1, -c*w_0*x_6*y_0 = 1, -z_2 = 1, -c*w_0*x_0*y_3 = 1, h*z_5 = 1, -z_5 = 1, c*q*w_3*y_1 = 1, y_6 = 1, c*q*w_6*y_0 = 1, -20*c*w_3*x_1*y_1 = 1, -3811631523823849533748331710598418339419473585110868329611049199088727487210518358698403600650604163616865365446917174642824731965169225010856 = 1, -6*c*w_2*x_0*y_2 = 1, -12*c*w_1*y_1 = 1, -815402048292665350361605655594682735739551400388698172957514912181357529452166069241702849357162652363127320586148655376530353341138636285780490937666199100656838591589818068069583942858123853509856470728047200160404503620011204187998022455209975923881995078157740323056620551393441867329220664704405892817577299688714602372 = 1, -5*c*w_0*x_1*y_4 = 1, c*q*w_5*y_1 = 1, x_3 = 1, c*q*w_2*y_2 = 1, -c*w_0*x_3*y_0 = 1, -w_5 = 1, beta*v_2*x_1 = 1, z_4 = 1, -6*c*q*w_2*y_2 = 1, -60*c*w_1*x_3*y_2 = 1, -3*c*w_1*x_0*y_2 = 1, -5*c*q*w_1*y_4 = 1, -7118526460094158566237096960490561709091070092828086062441195468912695928415522698155996215931724178282067881868879638064579962114818990405168966475255010638754980655611998225896646171876010907661249599240358565428108238285506433877468031680719737163294282252626338152191118292372037938699102061638984825319150824776777527388396550000914304469622110278259228804455923450618434869065040 = 1, -3*c*w_0*y_1 = 1, -beta*v_0*x_0 = 1, b*w_1 = 1, c*q*w_0*y_2 = 1, c*q*w_0*y_5 = 1, -10*c*w_3*x_0*y_2 = 1, y_4 = 1, -4*c*w_1*x_3*y_0 = 1, -18868275170634996523136416062498983968867351914014193414471986630184526794799440538880069229480783183353398157120903083393105059065818070513252870661353206952132444109638512536845548742781744068390199300 = 1, -35238011654326437371059416444119838700537121247060979719137014113665350279779162354060118065937189034372584167858558515279505959401029616894137584301857373169751419936952549627317562604131276306642254481919809910306853648744994676069811331658693417913129106079246372861100360171488943149337250394854532631205302046461561062085283811469131525517147743176523043922730426461827142223497192841862388006787187122123531501245381782281308422794350023432 = 1, -60*c*w_3*x_1*y_2 = 1, c*q*w_4*y_1 = 1, y_1 = 1, -z_1 = 1, a*y_5 = 1, -4*c*q*w_3*y_1 = 1, -10*c*w_2*x_3*y_0 = 1, c*q*w_1*y_2 = 1, -30*c*w_1*x_1*y_4 = 1, b*w_3 = 1, beta*v_0*x_1 = 1, -beta*v_0*x_1 = 1, h*z_3 = 1, u*v_1 = 1, c*q*w_5*y_0 = 1, w_5 = 1, -c*q*w_0*y_5 = 1, -beta*v_0 = 1, -6*c*w_0*y_2 = 1, -beta*x_0*v_1 = 1, beta*v_3*x_1 = 1, -10*c*w_0*x_3*y_2 = 1, a*y_1 = 1, c*q*w_1*y_5 = 1, d*x_5 = 1, beta*v_0*x_5 = 1, c*q*w_1*y_3 = 1, -3*c*w_1*y_0 = 1, c*q*w_4*y_2 = 1, y_5 = 1, -c*q*w_4*y_0 = 1, -w_7 = 1, z_5 = 1, b*w_5 = 1, v_1 = 1, h*z_1 = 1, z_aux = 1, -beta*v_0*x_3 = 1, -y_4 = 1, a*y_4 = 1, -c*w_0*x_0*y_4 = 1, x_1 = 1, -5*c*w_0*x_4*y_1 = 1, -5*beta*v_4*x_1 = 1, -10*c*q*w_2*y_3 = 1, -4*c*q*w_1*y_3 = 1, -w_2 = 1, beta*v_0 = 1, z_2 = 1, w_7 = 1, -4*beta*v_3*x_1 = 1, -30*c*w_2*x_1*y_2 = 1, -w_4 = 1, -15*c*w_4*x_0*y_2 = 1, -20*c*w_0*x_3*y_3 = 1, -2*beta*v_1*x_1 = 1, d*x_4 = 1, b*w_4 = 1, beta*v_1*x_1 = 1, a*y_3 = 1, c*q*w_0*y_0 = 1, -15*c*w_0*y_4 = 1, -12*c*w_2*x_1*y_1 = 1, h*z_4 = 1, -20*c*w_3*x_3*y_0 = 1, -c*q*w_3*y_0 = 1, -y_3 = 1, -4*c*w_0*x_1*y_3 = 1, -w_3 = 1, -6*c*w_1*x_5*y_0 = 1, d*x_1 = 1, -10*c*q*w_3*y_2 = 1, -c*w_1*x_0*y_0 = 1, beta*v_1 = 1, -beta*v_2*x_0 = 1, -30*c*w_4*x_1*y_1 = 1, v_3 = 1, -lm = 1, -6*c*w_1*x_1*y_1 = 1, -30*c*w_1*x_4*y_1 = 1, -c*q*w_0*y_4 = 1, -y_2 = 1, -10*c*w_2*x_0*y_3 = 1, -c*q*w_0*y_0 = 1, c*q*w_1*y_4 = 1, c*q*w_1*y_1 = 1, -beta*v_5*x_0 = 1, w_2 = 1, -3*c*q*w_1*y_2 = 1, h*z_0 = 1, -c*w_0*x_0*y_2 = 1, -5*c*q*w_4*y_1 = 1, -15*c*w_0*x_4*y_2 = 1, -c*q*w_1*y_0 = 1, -3*c*q*w_2*y_1 = 1, -c*q*w_2*y_0 = 1, w_3 = 1, z_3 = 1, -c*w_0*x_5*y_0 = 1, beta*v_4*x_1 = 1, y_2 = 1, -436608920299180809179630557104527219900824252088990974277221602192823366649534158 = 1, -2*c*q*w_1*y_1 = 1, beta*v_4*x_0 = 1, w_6 = 1, -4*c*w_0*x_3*y_1 = 1, beta*v_2*x_0 = 1, -20*c*w_1*x_1*y_3 = 1, x_4 = 1, -c*w_0*y_0 = 1, -c*q*w_5*y_0 = 1, -5*beta*v_1*x_4 = 1, x_5 = 1, -2*c*w_0*x_1*y_1 = 1, -c*w_3*x_0*y_0 = 1, c*q*w_2*y_1 = 1, -w_0 = 1, beta*v_0*x_3 = 1, -c*q*w_0*y_3 = 1, d*x_0 = 1, beta*v_0*x_4 = 1, -90*c*w_2*y_2 = 1, -5*c*w_4*x_1*y_0 = 1, -4*c*w_3*x_0*y_1 = 1, -10*c*w_3*y_0 = 1, v_2 = 1, -z_0 = 1, v_3*beta = 1, -5*c*w_1*x_4*y_0 = 1, b*w_0 = 1, -60*c*w_3*y_1 = 1, -4*c*w_3*x_1*y_0 = 1, y_3 = 1, c*q*w_1*y_0 = 1, -z_4 = 1, -10*beta*v_2*x_3 = 1, -c*w_0*x_0*y_6 = 1, -5*c*w_1*x_0*y_4 = 1, c*q*w_2*y_0 = 1, c*q*w_2*y_3 = 1, -6*c*w_5*x_1*y_0 = 1, -z_3 = 1, -c*w_0*x_0*y_1 = 1, -60*c*w_1*y_3 = 1, -c*w_5*x_0*y_0 = 1, -60*c*w_2*x_1*y_3 = 1, -3*c*w_2*x_1*y_0 = 1, -c*w_0*x_0*y_0 = 1, w_4 = 1, beta*v_0*x_0 = 1, -c*q*w_0*y_2 = 1, -beta*v_4*x_0 = 1, c*q*w_4*y_0 = 1, d = 1, -2*c*w_1*x_0*y_1 = 1, -6*c*w_0*x_1*y_5 = 1, -c*w_6*x_0*y_0 = 1, -w_1 = 1, -y_1 = 1, -c*w_4*x_0*y_0 = 1, -c*w_0*x_1*y_0 = 1, u*v_0 = 1, beta*x_0*v_1 = 1, -6*c*w_5*x_0*y_1 = 1, c*q*w_2*y_4 = 1, -20*c*w_1*x_3*y_1 = 1, beta*v_5*x_0 = 1, -30*c*w_2*y_1 = 1, -w_6 = 1, -3*c*w_2*x_0*y_1 = 1, -1 = 1, -15*c*w_2*x_0*y_4 = 1]
# 273, -5.537035820
# 1
# [k = 5, x_2 = 25]
# [k = [2, 2, 2, 2, 2], x_2 = [4, 1, 4, 4, 3, 2, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 3]]