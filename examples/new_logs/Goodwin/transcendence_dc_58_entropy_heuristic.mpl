infolevel[Groebner]:=10;
Et_hat := [1287373283726616403-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-5583785807904221688091397891714113839806455028860598639/8060197978772174256, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+15416636963751088380930499081016289060802167651805441435384754399301432391226793640237465102564038169287629533/41355002868899505083038687868856677219937118939875501312, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-21282397234476190126400806538782783188776528615903695501226454706281166458200676721343841310902708384199244265032348825027696517738458039381581403058967056778745279/106091454998429823638989533444445495369863949730374032168814906111611176993718892730096474112, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4-4469985352343876984975589550950381500743244706620398923026095462480295612504983878704616400662716996564845124614927204296913567571677514071568550859586290515445383611766566997736775232270068856568643616798660576883829143771/1632991808903730385551738165705220035758562572185240937550751787591255702092176399841856619978431092904794406510016173951866601472, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-447827748406383570500457188342785108822800953087406059147798056953635595373326217600663969332801736321231879306163556807533925522428413190448299694224057088640779422719636293110986195837187669606395289558615742646429952102056335934021361330566002683198711208371368845477080929273840435799982196165/2094625377688008797338020421783022638150218340872247744155592290747776514513403117406143514839464539611649430655602949289410849879324340792527231407366196471282139136, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-111571974804474537117579602621858361248743616707603383712669455725452699721360408060837108890525237730529256772080148234194505898581996137759917582404628435667898764414818191578503543805089552924182209253726235450969608024699292563662430199959569122325789108377521865862113141361954962279558726702278379794622158118025465715231869694587578918327792667247761500891748939687/10747036081705874024974088808353772466480238467612260343013482604054408777537595181412098006731546826063070787430827111374810097645479879308123269255417876112197498820655792932481051609725036188319875072, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+2660333111930947804023207916152089022346527347354882617151895593250503815494169706324243795166956101569450810332254624431693480523946547880484052308975352031443063370334505286576725173607766539460722387894724100340548650462885842606703187493292442433129782054308959099222483391421544268861155341211803216570184757464271241507484582820234214919918621100097199784819134888291460104435764496910536776227530708612924283245835893503372988068448491079/3063363601693895480188293553981478131498967052932229955830161385202276860829362948970219852741324536272277473219233198236025561297467903142954499457497798059642802650420592680438562566112775135026277734002795184166408732192903241775710208, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+2599322640150031028340022654109759974562318832068529848196798072683807490393029105824189908082100492987100464018185094734760480478072823107779775666387758812177215988418118897360784517286946557272673911112307002591079494944335654164934768669635959478692612498526349914246329873104066802763063919900593867455917514880456603736766654951792063206748841096197337564374772920140060191132778901284559800474568221880385575036414341355724856523379850483764060961427135489616620868648269394500128010301949884844431671487794441513/6736031486545531216372390352261043380211735200420488134522827602058891376765346451103576605407081173441442783959523579311618493187334845842206635670944703472896143139026103895824412362064197391622424453405852649531297974965164911299196794436520870000875422740596221294936064, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+20597883710168462498358034994615801390142538369913572354930776076445482129128354298533739145689423929584491031665750043046428088795444742252315143233849983920845737899023072664265703832686283275820689593418420393941242133588413798867729105465959464000160845159953132469183181919354078500554679926811645557545297755888173208099476008175633423946211737058817564606934082973551459947075714083813246988517736890113687399564260771917926154511300943904363236536071983107799261262125397460002289552485177065137675794113197367355735140971302546236757051127405367922110051800523117539629572162294205985/270007908561092989724037292099124150403045937675547909040449882038897227936567480711128844094768266293537676844832158296316960204317475382558489444119734973022347426600641028860606725724816515739383882618036370671922577144950257167051367681031787204248854104062713244845396473770948287675249666455354554187776, 3058419358542345187309416401863070725-x2_1, 1319033907165593891306098883571105572312898306308776337574280318820898397266726373297405122247/8060197978772174256-x3_3, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_1, x3_3] 66
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, x3_2*x4_3 = 1, b*x1_6*x4_0 = 1, x4_1 = 1, x3_7 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, delta*sgm*x3_2*x4_3 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, x3_2 = 1, x3_6*x4_1 = 1, x2_3 = 1, b*x1_2*x4_2 = 1, -gama*x2_4 = 1, x3_2*x4_5 = 1, delta*x3_4 = 1, x3_2*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, -3*gama*x4_2*sgm = 1, delta*sgm*x3_6*x4_0 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -gama*sgm*x2_3*x4_0 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, -x1_7 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x2_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, x2_4 = 1, -2*gama*x4_1*sgm = 1, -gama*x2_3 = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x3_4*x4_2 = 1, x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -gama*sgm*x2_2*x4_0 = 1, -15*gama*sgm*x2_4*x4_2 = 1, b*c*x1_7 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, -7*gama*x4_6*sgm = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, x2_5 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, delta*x3_2 = 1, -gama*sgm*x2_0*x4_0 = 1, -20*gama*sgm*x2_3*x4_3 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, x3_7*x4_1 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, -6*gama*x4_5*sgm = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, -gama*sgm*x2_4*x4_0 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, x4_3 = 1, x4_5 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, -4*gama*x4_3*sgm = 1, x1_2*x4_2 = 1, x3_5 = 1, -5*gama*x4_4*sgm = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, x4_4*delta*sgm = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, x4_1*delta*sgm = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, x1_3*x4_2 = 1, beta = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, delta*sgm*x3_2*x4_2 = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x1_3*x4_1 = 1, -gama*x2_6 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, -x1_3 = 1, x4_2*delta*sgm = 1, -10*gama*sgm*x2_3*x4_2 = 1, x3_2*x4_6 = 1, -x1_5 = 1, -gama*x2_2 = 1, delta*sgm*x3_1*x4_1 = 1, x4_3*delta*sgm = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, -21*gama*sgm*x2_2*x4_5 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_4*x4_3 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x4_4 = 1, x3_0*x4_7 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, x3_2*x4_1 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, delta*sgm*x3_2*x4_5 = 1, -gama*x2_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, -gama*sgm*x2_0*x4_5 = 1, delta = 1, x1_1*x4_5 = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1, -gama*sgm*x2_5*x4_0 = 1]
# 282, -5.577999884
# 2
# [x2_1 = 10, x3_3 = 12]
# [x2_1 = [4, 1, 4, 2, 2, 4, 4, 4, 4, 4], x3_3 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2]]
# [x2_1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], x3_3 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]