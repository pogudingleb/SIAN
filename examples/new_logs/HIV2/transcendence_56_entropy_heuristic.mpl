infolevel[Groebner]:=10;
Et_hat := [228683851353891543799-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 153027932215937267319-z_0, -c*q*w_0*y_0+h*z_0+z_1, 196838919744990471510-k_0, k_1, -w_1+2073490763801341587654924927714151188285120656004197679375726264792385575606830684, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 1182730557958186477700648237715127242065030074736165370908966892458125006753030917-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+18800470265458943505684493615125023584660541883741503028853203007224251324035341198271471566610457635190446055191451838313033606426696404535332, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, 10723891842265732813965393118801197224376607684150113571104612264666609695850993855324061722698739560665072110023481910649637991494192594054219-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+170465038172829879811692369379900393321430974190023531469906546852893115520402656760146608355861249964261864093590931656486546003543948119873298053088355609108575304300692173535039376391248987919071067648, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 97234197147275564978170405845436266319796409802273127802495226450037269426372331903235154642307528460161502077732585724713740532852305187064973803822610840691280963831901010823354938331264463281434735181-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+1545617148346102667907273876045169430622749092863526665993619969271221462046166730581559735954866601194762779951523061882345739977746925420123455133766077965370165011359456132779435851909387637654812039960313650706258252024818950155916667442748560991317129606839736, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, 881628538774754828971256464738022015918260465282827021189700514309425277023458120112385280124779651329489686994403401759886851510721328151881796337037855609444557300746544600548275222616413208530846931200935626920097850998488235137816230924486242394320066163861507-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5+14014207223181240196906775989984133945673675270378366229739428980214108231825344781631885866226459317166890434666169074937302531578810872052074608917193646416556062243181959676043192409640939268059373463547446606505450192373225676031540440369213133038167461151118103782640552447820116205131696816089753826653329049704710043680, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 7993781027520808978439256352765536731742617088458246017581012936116214676267119514248065414032674688503458158256669377880539701224705012899157052946932205189618724558106789357804809011143906773358009895285836287875606584796517704179698071212019267632888099567001879391413615564386854066580259608453614503973076712786136688573-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+127067692218879794545602056441213042175528634735690961208446395548435719151474654574898372674353284847997156398840392391678706464570233186621090692005504948865666099824596322471978542062978142504589744278307400969273794491146918893121869922902931712853277048300881319540782093315015386962977389401690442959515411012045729982442799475216622243262086517776433395476192688533852182764559008, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7+1152130701986777114556058756926317036908422693529438305337323759053086410098077522026480712522946464224352999188496180725058043207608856197915283359312468982475706112437858841541600912904813571792305185370246361035255324939898969701214128847898694443748994039170563270326127997145437435510117742833577932723281793300053759482760147883781142445008285408782773617549709821415743946598969496982541928981195970303510735993741663480210355585730985227136, 72480111867473740870029720314695324140465027418060733531722249426542248508605741342570071184348578558682074594139847643203408852904438687390774321800293648238711378552806124408327514050503525573229645470545926094681796564576679840974214145493258772019574621721666395494040301610285989287636634055180826574197308552310606402805495219750668166651256105019382592166240804555731406296711587-z_6, -k_1, -k_2, -k_3, -k_4, 277147729972003337785-x_0, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, x_0] 160
# 266, [1 = 8, -15*c*w_2*x_4*y_0 = 1, c*q*w_3*y_2 = 1, -3*c*w_0*x_2*y_1 = 1, -164818727526407584224582241015472061561268250705468742805250854340282119696147721946275673955778194068305854889797040849026418884917373978479502188708810660198213701132365888216757381609938143300129045980020666701804658690708597812560949912633477767356876233539525765096819262942688165855214619303333905173103926951655537747511517048619957434971519104080160061107503551598056956786627727763745356294986395603386928265724434770341498707000415 = 1, -5*c*w_4*y_1 = 1, -30*c*w_1*x_1*y_4 = 1, v_3*beta = 1, c*q*w_2*y_1 = 1, beta*v_1 = 1, w_7 = 1, -z_3 = 1, -6*c*w_5*x_1*y_0 = 1, a*y_1 = 1, -c*w_0*y_2 = 1, z_aux = 1, -beta*v_0*x_5 = 1, c*q*w_3*y_3 = 1, -278915016615161777466941150036983400437558304405153015500624043150453025912641951665247663197050419117810744691810756103357604551176032186838 = 1, -beta*v_0*x_4 = 1, -10*c*w_3*y_2 = 1, v_4*beta = 1, -20*c*w_1*x_3*y_1 = 1, v_2*beta = 1, -v_3*beta = 1, w_1 = 1, -z_6 = 1, -v_5*beta = 1, -4*c*w_3*y_1 = 1, -2*c*w_1*x_1*y_0 = 1, -y_1 = 1, -2*c*q*w_1*y_1 = 1, a*y_5 = 1, -5*c*w_1*x_4*y_0 = 1, -z_0 = 1, beta*v_0*x_5 = 1, c*q*w_3*y_1 = 1, beta*v_0*x_4 = 1, beta*v_4*x_1 = 1, -z_4 = 1, x_5 = 1, beta*v_3*x_1 = 1, -c*w_5*y_0 = 1, h*z_3 = 1, a*y_2 = 1, -10*c*w_3*x_2*y_0 = 1, w_3 = 1, -c*q*w_4*y_0 = 1, y_4 = 1, d*x_3 = 1, -c*q*w_0*y_1 = 1, -4*beta*v_1*x_3 = 1, -y_3 = 1, c*q*w_3*y_0 = 1, v_5*beta = 1, -60*c*w_1*x_2*y_3 = 1, z_1 = 1, -c*w_4*y_0 = 1, z_4 = 1, -c*q*w_0*y_3 = 1, -60*c*w_2*x_3*y_1 = 1, beta*v_0*x_3 = 1, -c*w_0*x_3*y_0 = 1, -c*q*w_0*y_2 = 1, -60*c*w_2*x_1*y_3 = 1, -6*c*w_1*x_1*y_1 = 1, -3*beta*v_2*x_1 = 1, beta*v_1*x_1 = 1, x_6 = 1, -6*beta*v_2*x_2 = 1, -5*c*w_1*y_4 = 1, -c*q*w_1*y_0 = 1, -w_6 = 1, -30*c*w_1*x_4*y_1 = 1, -3*c*w_2*x_1*y_0 = 1, -10*beta*v_2*x_3 = 1, h*z_5 = 1, -z_2 = 1, u*v_4 = 1, beta*v_0*x_1 = 1, -6*c*w_0*x_1*y_5 = 1, d*x_2 = 1, c*q*w_6*y_0 = 1, -10*beta*v_3*x_2 = 1, w_4 = 1, c*q*w_1*y_2 = 1, -c*w_0*y_3 = 1, -15*c*w_4*x_2*y_0 = 1, c*q*w_2*y_4 = 1, u*v_2 = 1, y_5 = 1, -4*c*w_3*x_1*y_0 = 1, -90*c*w_2*x_2*y_2 = 1, z_5 = 1, c*q*w_1*y_0 = 1, -60*c*w_1*x_3*y_2 = 1, beta*v_2*x_1 = 1, -c*w_0*x_1*y_0 = 1, -15*c*w_2*y_4 = 1, -c*w_0*x_2*y_0 = 1, c*q*w_0*y_0 = 1, -2*c*w_0*x_1*y_1 = 1, d = 1, -5*beta*v_1*x_4 = 1, -c*q*w_5*y_0 = 1, -w_3 = 1, -w_0 = 1, -20*c*w_1*x_1*y_3 = 1, beta*v_0 = 1, -c*w_0*y_6 = 1, c*q*w_4*y_0 = 1, -6*c*w_1*x_5*y_0 = 1, -beta*v_0 = 1, beta*v_1*x_3 = 1, -5*c*w_4*x_1*y_0 = 1, beta*v_3*x_2 = 1, c*q*w_0*y_4 = 1, -3*c*w_1*x_2*y_0 = 1, -4*beta*v_3*x_1 = 1, -w_5 = 1, -6*c*w_1*y_5 = 1, a*y_4 = 1, z_6 = 1, u*v_0 = 1, h*z_2 = 1, -30*c*w_1*x_2*y_2 = 1, -10*c*w_2*y_3 = 1, -w_2 = 1, beta*v_1*x_2 = 1, -3*c*w_1*y_2 = 1, c*q*w_1*y_1 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_1*y_4 = 1, -15*c*w_0*x_2*y_4 = 1, c*q*w_0*y_3 = 1, -c*q*w_0*y_4 = 1, -10*c*q*w_3*y_2 = 1, -6*c*w_0*x_2*y_2 = 1, -4*c*w_1*y_3 = 1, -c*w_0*x_5*y_0 = 1, c*q*w_4*y_1 = 1, -y_0 = 1, -6*c*w_5*y_1 = 1, x_4 = 1, beta*v_0*x_2 = 1, c*q*w_5*y_1 = 1, -10*c*w_2*x_3*y_0 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, -15*c*w_4*y_2 = 1, u*v_1 = 1, c*q*w_1*y_5 = 1, y_6 = 1, -lm = 1, -3*c*w_2*y_1 = 1, h*z_1 = 1, -4*c*q*w_1*y_3 = 1, -6*c*w_2*x_2*y_0 = 1, -y_0*w_0*c = 1, -20*c*w_0*x_3*y_3 = 1, -2*beta*v_1*x_1 = 1, -5*c*w_0*x_4*y_1 = 1, y_2 = 1, -200275591362100408122726486434978192082818297638257436524234088172441429289786695546364547781045103911678228724741854185416225592368885577976888230304127670336166398780833281003232908235933820396976353896937396568637104552052874068109358790706811413731355547934775775296348502446628354538461281317537642566573645375787349764862127053001630327821232612057741174355984636146724963318 = 1, -20*c*w_3*x_1*y_1 = 1, -10*c*q*w_2*y_3 = 1, d*x_4 = 1, -4*c*w_0*x_3*y_1 = 1, -5*beta*v_4*x_1 = 1, a*y_3 = 1, beta*v_2*x_2 = 1, c*q*w_1*y_3 = 1, x_1 = 1, w_5 = 1, -w_4 = 1, beta*v_2*x_3 = 1, -194503980712306719331673225396981017649589001127842111128364091688806229397220549852914261249712612573911644372217205957552089190787944608959325475955080188597633654975933405822389120551873648823906260156158992646533360385166232031578306153955993131239259904604743847353603901133386784210582875258498557781532618530608495 = 1, -beta*v_0*x_3 = 1, -12*c*w_2*x_1*y_1 = 1, -c*q*w_2*y_0 = 1, c*q*w_0*y_2 = 1, -15*c*w_0*x_4*y_2 = 1, c*q*w_0*y_5 = 1, -beta*v_0*x_2 = 1, -3*beta*v_1*x_2 = 1, z_3 = 1, -20*c*w_3*x_3*y_0 = 1, -6*c*w_2*y_2 = 1, -c*w_3*y_0 = 1, -c*w_0*y_5 = 1, -20*c*w_3*y_3 = 1, b*w_2 = 1, -10*c*w_0*x_3*y_2 = 1, b*w_0 = 1, -beta*v_1 = 1, -60*c*w_3*x_2*y_1 = 1, b*w_3 = 1, d*x_1 = 1, -6*c*w_0*x_5*y_1 = 1, z_2 = 1, -4*c*q*w_3*y_1 = 1, y_1 = 1, -30*c*w_2*x_1*y_2 = 1, -2*c*w_1*y_1 = 1, -w_7 = 1, -c*w_0*y_1 = 1, v_5 = 1, -4*c*w_1*x_3*y_0 = 1, x_2 = 1, -270877148049476526865509287601522425205673947385712555130824215342753506329726095 = 1, c*q*w_5*y_0 = 1, c*q*w_4*y_2 = 1, c*q*w_2*y_3 = 1, -5*c*q*w_1*y_4 = 1, -c*w_0*x_4*y_0 = 1, -c*q*w_0*y_0 = 1, b*w_6 = 1, -v_2*beta = 1, c*q*w_2*y_0 = 1, c*q*w_0*y_6 = 1, -229535800213430833046071428647745561111868539496626746574080018173889401798672933180034287262201265441088228325350952475773216169046077712398889968292300765668791503946947052542571478451556870450681215 = 1, u*v_3 = 1, -z_1 = 1, -3*c*q*w_1*y_2 = 1, -z_5 = 1, w_6 = 1, v_4 = 1, -10*c*w_0*x_2*y_3 = 1, a*y_0 = 1, -5*c*w_0*x_1*y_4 = 1, h*z_0 = 1, v_1 = 1, -236346926978904436719568454748778602012358613691752932997434020712583894484443859566116894969380093568709675867745002027093398462799396399447923208012074320225663661281767874529591683900145755749483687450867533603163647030756744060339181405835628716291944151078 = 1, h*z_4 = 1, -12*c*w_1*x_1*y_2 = 1, -c*w_2*y_0 = 1, -w_1 = 1, -30*c*w_2*x_2*y_1 = 1, w_2 = 1, -c*w_0*x_6*y_0 = 1, b*w_1 = 1, -c*w_6*y_0 = 1, c*q*w_2*y_2 = 1, -c*q*w_3*y_0 = 1, b*w_4 = 1, beta*v_1*x_4 = 1, -60*c*w_3*x_1*y_2 = 1, -c*w_0*y_4 = 1, -12*c*w_1*x_2*y_1 = 1, -c*w_1*y_0 = 1, d*x_5 = 1, v_3 = 1, -3*c*w_0*x_1*y_2 = 1, -y_4 = 1, v_2 = 1, b*w_5 = 1, -3*c*q*w_2*y_1 = 1, y_3 = 1, c*q*w_0*y_1 = 1, -6*c*q*w_2*y_2 = 1, -v_4*beta = 1, -beta*v_0*x_1 = 1, -y_2 = 1, x_3 = 1, -1 = 1, -c*q*w_0*y_5 = 1]
# 273, -5.548535779
# 1
# [k = 5, x_0 = 41]
# [k = [2, 2, 2, 2, 2], x_0 = [4, 4, 4, 3, 2, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3]]