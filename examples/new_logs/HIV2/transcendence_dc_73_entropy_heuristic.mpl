infolevel[Groebner]:=10;
Et_hat := [115143198838659752517-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 484017552638891914596-z_0, -c*q*w_0*y_0+h*z_0+z_1, 211037073962289384411-lm_0, lm_1, -w_1-87502796468490656470451459488323539461506115826324625804093955539842970791477118, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm_0+x_1, -beta*v_0*x_0+a*y_0+y_1, 1781710859339938981194595170957937740052006509782265115252094202392523142313765202-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+66497539299171550182449674703488504631944733220844914513826739525943668347177465150106824412844395596850243261482544958207169822946992217922, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1-lm_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k*y_0+u*v_0+v_1, -1354006873727540517157257575916001280644182614306144000784197090793749310102720399474998435640500678243739116829923562597540971066306941635726-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3-50534644734893452416293921332908906507531098608909745343471408104645920309207635510730232091216019254111834336256364739295018438024892154888091501486058178314277895808583083096451889686398930123651988, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2-lm_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, lm_2, -k*y_1+u*v_1+v_2, 1028974260605120657363868585883073701638849340196940816600638678373711354410751295253572385464788783661982229713827744281913967678935998087817102483977512442410004888576152250498228468835729622812641738-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+38403681480492166885188232442662573184044269989982262412272934743248478505474895175877698259225742055982220427466696680875191671818634233377808717020634123609829536891305586584197391827958854627625458951770897376961927164425239551959407133119502109194121017402, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4-lm_3, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, lm_3, -k*y_2+u*v_2+v_3, -781966509573945450180197202709686476795482945614105975807317685819373151194578017839104129404536595513447491093266095939472542627182952480253640289288827854458026661110253826492154076903433794634997385018717031109957190344902390454071591983879184913421945057994-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5-29184785190282323593602664316823801256187535729716136767115411977060627309684441824404320573684688054092911303251830862855686514984237258238194369638417649888128621474141740178449355956845342835329036285006401333710888654138915314014696705643861730751482401553611672185934124911455539072584179740283759596470098265762408, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5-lm_4, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, lm_4, -k*y_3+u*v_3+v_4, 594253564453268454870198861557722303620785031344116242442340792946077180036625238666961629059313364477638131792778294844899917680767011305235993744251569060461530057388298067133236844074829681482058591426330356137060766716796525015589933057177625676823496717554294334183063116277839965136949137396219556218626945708531622-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+22178907171584184226745615817920416796288090111107598261828467601034873042391473382129852931348845262698222567828514942349202402219746476064778057416359684162210622954360065765175180296412565096264379254024402561896553253833382789582457264315912698365929826380685147548860377749173836585340196702295160883087234341083341738878226883349756993545537875801735394713614618853934784082, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6-lm_5, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, lm_5, -k*y_4+u*v_4+v_5, -w_7-16854807055065730579366655645715029110677448718220330802574906206634992330020170805054717777725229440904507728944183663274704645195124368412139211352400075547853436961211291888485548891187907106703729154420346147888747303047736851978395878342095415346844048948776443114511505939400712977944867215542859167894740621422121794019747102423032250965040211647834686009014319775625554678180712178233879230012203030745464065446391156497963439874128, -451601564186964715819323793136017659235763679396280219990101984363904948039916912898844762987777095504691610427996191942009001466272420616890084228532638705044245087452434980813639894490749169707346372904924884419339791951318617749688989379948386278637196458685077614412222913079836618420878919931557954466948164929065350805882541499496381028492046790172965926194740947412978553686-z_6, -lm_1, -lm_2, -lm_3, -lm_4, -lm_5, 1470497154606070473546186699306050264996593223928670156111892876642054895889535661-v_2, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, lm_5, z_5, y_5, x_5, w_5, v_5, lm_4, z_4, y_4, x_4, w_4, v_4, lm_3, z_3, y_3, x_3, w_3, v_3, lm_2, z_2, y_2, x_2, w_2, v_2, lm_1, z_1, y_1, x_1, w_1, v_1, lm_0, z_0, y_0, x_0, w_0, v_0, z_aux, w_aux, a, b, beta, c, d, h, k, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [lm, v_2] 28
# 257, [1 = 16, -1 = 2, v_4 = 1, a*y_0 = 1, b*w_1 = 1, -10*c*w_3*x_0*y_2 = 1, -w_0 = 1, -6*c*w_2*x_2*y_0 = 1, -2*c*w_1*x_0*y_1 = 1, -4*c*q*w_1*y_3 = 1, -15*c*w_4*x_2*y_0 = 1, -c*w_1*x_0*y_0 = 1, c*q*w_2*y_3 = 1, c*q*w_0*y_2 = 1, -5*c*w_1*x_0*y_4 = 1, beta*v_0*x_5 = 1, -c*w_0*x_0*y_4 = 1, -beta*v_0*x_3 = 1, -c*w_0*x_5*y_0 = 1, c*q*w_1*y_2 = 1, u*v_4 = 1, beta*v_3*x_2 = 1, -c*q*w_5*y_0 = 1, w_1 = 1, -k*y_3 = 1, beta*v_3*x_1 = 1, -15*c*w_0*x_4*y_2 = 1, -3*c*w_1*x_2*y_0 = 1, -k*y_4 = 1, -c*w_0*x_1*y_0 = 1, -4*c*w_3*x_0*y_1 = 1, -2*c*w_1*x_1*y_0 = 1, -10*c*q*w_3*y_2 = 1, c*q*w_3*y_1 = 1, h*z_2 = 1, -6*c*w_0*x_5*y_1 = 1, -3*c*w_2*x_0*y_1 = 1, -6*c*w_0*x_2*y_2 = 1, b*w_6 = 1, x_3*beta = 1, -12*c*w_2*x_1*y_1 = 1, y_2 = 1, -c*w_0*x_2*y_0 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_5*y_1 = 1, -c*w_3*x_0*y_0 = 1, -3*beta*v_1*x_2 = 1, c*q*w_2*y_4 = 1, -z_5 = 1, -15*c*w_4*x_0*y_2 = 1, -3*c*w_2*x_1*y_0 = 1, -beta*v_5*x_0 = 1, b*w_2 = 1, -c*w_0*x_0*y_1 = 1, v_1 = 1, -z_3 = 1, -w_6 = 1, -5*c*w_0*x_4*y_1 = 1, -c*q*w_0*y_4 = 1, beta*v_1*x_2 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, beta*v_0*x_1 = 1, -20*c*w_3*x_1*y_1 = 1, -c*w_0*x_0*y_0 = 1, -c*w_0*x_6*y_0 = 1, b*w_3 = 1, z_2 = 1, y_3 = 1, -z_0 = 1, c*q*w_5*y_0 = 1, x_6 = 1, -10*c*q*w_2*y_3 = 1, d*x_2 = 1, z_5 = 1, -3*beta*x_1 = 1, -c*q*w_0*y_2 = 1, -c*w_0*x_3*y_0 = 1, -c*q*w_0*y_3 = 1, beta*v_1*x_1 = 1, -10*c*w_0*x_3*y_2 = 1, u*v_1 = 1, beta*v_1*x_4 = 1, -30*c*w_2*x_1*y_2 = 1, w_7 = 1, a*y_3 = 1, beta*x_0*v_1 = 1, -2*c*q*w_1*y_1 = 1, -w_2 = 1, z_4 = 1, -beta*x_0 = 1, -20*c*w_3*x_3*y_0 = 1, c*q*w_2*y_1 = 1, -c*q*w_4*y_0 = 1, -6*c*w_1*x_0*y_5 = 1, -w_5 = 1, -6*c*w_2*x_0*y_2 = 1, y_1 = 1, -6*c*q*w_2*y_2 = 1, c*q*w_2*y_2 = 1, -w_7 = 1, -10*c*w_0*x_2*y_3 = 1, -3*c*w_1*x_0*y_2 = 1, -10*c*w_2*x_3*y_0 = 1, -12*c*w_1*x_2*y_1 = 1, d*x_4 = 1, -60*c*w_3*x_2*y_1 = 1, c*q*w_6*y_0 = 1, -c*q*w_0*y_0 = 1, c*q*w_1*y_5 = 1, -c*q*w_0*y_5 = 1, -c*q*w_1*y_0 = 1, a*y_1 = 1, -4*c*w_1*x_3*y_0 = 1, -c*q*w_0*y_1 = 1, c*q*w_0*y_1 = 1, -3*c*q*w_1*y_2 = 1, beta*v_0*x_2 = 1, -30*c*w_1*x_1*y_4 = 1, x_3 = 1, -5*c*w_0*x_1*y_4 = 1, -w_3 = 1, -5*c*q*w_1*y_4 = 1, -c*w_2*x_0*y_0 = 1, -4*c*q*w_3*y_1 = 1, -30*c*w_2*x_2*y_1 = 1, h*z_4 = 1, -c*q*w_2*y_0 = 1, -6*c*w_1*x_1*y_1 = 1, -20*c*w_1*x_3*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_4 = 1, c*q*w_0*y_4 = 1, -3*c*w_0*x_1*y_2 = 1, -20*c*w_3*x_0*y_3 = 1, -5*beta*v_1*x_4 = 1, b*w_4 = 1, -z_6 = 1, c*q*w_4*y_2 = 1, -4*c*w_1*x_0*y_3 = 1, v_5 = 1, h*z_0 = 1, beta*v_0*x_3 = 1, -beta*v_0*x_1 = 1, c*q*w_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -beta*v_0*x_2 = 1, -c*w_0*x_0*y_6 = 1, c*q*w_4*y_0 = 1, x_2*beta = 1, -60*c*w_1*x_2*y_3 = 1, c*q*w_3*y_0 = 1, -12*c*w_1*x_1*y_2 = 1, beta*v_5*x_0 = 1, c*q*w_0*y_3 = 1, b*w_5 = 1, -c*w_0*x_0*y_3 = 1, x_2 = 1, beta*v_1*x_3 = 1, -5*c*w_1*x_4*y_0 = 1, -c*q*w_3*y_0 = 1, z_1 = 1, -5*c*w_4*x_1*y_0 = 1, -60*c*w_3*x_1*y_2 = 1, beta*v_0*x_4 = 1, -w_4 = 1, -w_1 = 1, c*q*w_3*y_2 = 1, -beta*x_0*v_1 = 1, beta*v_0*x_0 = 1, -30*c*w_1*x_4*y_1 = 1, x_1 = 1, beta*v_4*x_0 = 1, -10*c*w_3*x_2*y_0 = 1, -15*c*w_2*x_4*y_0 = 1, c*q*w_1*y_1 = 1, c*q*w_0*y_5 = 1, -4*beta*v_1*x_3 = 1, -3*c*w_0*x_2*y_1 = 1, -10*x_3*beta = 1, -5*c*w_4*x_0*y_1 = 1, -z_1 = 1, -z_4 = 1, beta*x_1 = 1, y_6 = 1, h*z_5 = 1, -10*c*w_2*x_0*y_3 = 1, c*q*w_0*y_0 = 1, x_4 = 1, beta*v_3*x_0 = 1, -beta*v_0*x_5 = 1, -60*c*w_1*x_3*y_2 = 1, w_5 = 1, w_2 = 1, d*x_1 = 1, -6*c*w_5*x_0*y_1 = 1, -6*c*w_1*x_5*y_0 = 1, -beta*v_3*x_0 = 1, c*q*w_1*y_3 = 1, h*z_1 = 1, h*z_3 = 1, w_3 = 1, c*q*w_3*y_3 = 1, d*x_5 = 1, v_3 = 1, -beta*v_0*x_0 = 1, -2*c*w_0*x_1*y_1 = 1, -20*c*w_0*x_3*y_3 = 1, x_5 = 1, b*w_0 = 1, -k*y_0 = 1, -4*beta*v_3*x_1 = 1, -4*c*w_0*x_3*y_1 = 1, -15*c*w_0*x_2*y_4 = 1, beta*v_4*x_1 = 1, -c*w_4*x_0*y_0 = 1, -5*beta*v_4*x_1 = 1, y_5 = 1, a*y_2 = 1, u = 1, -20*c*w_1*x_1*y_3 = 1, y_4 = 1, w_6 = 1, -60*c*w_2*x_3*y_1 = 1, -k*y_1 = 1, u*v_0 = 1, -c*w_0*x_0*y_2 = 1, a*y_5 = 1, -30*c*w_1*x_2*y_2 = 1, w_4 = 1, a*y_4 = 1, c*q*w_4*y_1 = 1, -z_2 = 1, -60*c*w_2*x_1*y_3 = 1, -2*beta*v_1*x_1 = 1, -c*w_5*x_0*y_0 = 1, -4*c*w_3*x_1*y_0 = 1, -6*c*w_5*x_1*y_0 = 1, -beta*v_4*x_0 = 1, -c*w_0*x_4*y_0 = 1, -k*y_2 = 1, d*x_0 = 1, z_3 = 1, d*x_3 = 1, u*v_3 = 1, -90*c*w_2*x_2*y_2 = 1, z_6 = 1, c*q*w_1*y_0 = 1, -6*x_2*beta = 1, -c*w_6*x_0*y_0 = 1, c*q*w_1*y_4 = 1, -6*c*w_0*x_1*y_5 = 1, c*q*w_2*y_0 = 1, -10*beta*v_3*x_2 = 1, beta*x_0 = 1, z_aux = 1, -15*c*w_2*x_0*y_4 = 1]
# 273, -5.441897751
# 2
# [lm = 1, v_2 = 10]
# [lm = [1], v_2 = [3, 3, 1, 3, 3, 2, 3, 3, 3, 3]]
# [lm = [1], v_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]