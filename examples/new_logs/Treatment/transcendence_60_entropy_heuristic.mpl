infolevel[Groebner]:=10;
Et_hat := [50583945842865-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 55327731963622-b_0, b_1, -Tr_1-1009829219590573693145316945, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b_0*d+In_0*N*a+In_0*N*g-In_0*S_0*b_0+In_1*N, -Tr_2+31548505119965332474031669495879934827779755447363317396930288390595/2753953988527, -In_2*g+Tr_2*nu+Tr_3, ((-Tr_1*d-In_1)*b_0-b_1*(Tr_0*d+In_0))*S_0-S_1*(Tr_0*d+In_0)*b_0+N*((a+g)*In_1+In_2), S_0*Tr_0*b_0*d+In_0*S_0*b_0+N*S_1, -Tr_3-498558625388671447592988021945172703713283448420977797535304439567051595877491544094745092349621143333671285/7584262570923771647629729, -In_3*g+Tr_3*nu+Tr_4, ((-Tr_2*d-In_2)*b_0+(-Tr_0*b_2-2*Tr_1*b_1)*d-In_0*b_2-2*In_1*b_1)*S_0+((-2*S_1*Tr_1-S_2*Tr_0)*d-S_2*In_0-2*S_1*In_1)*b_0-2*S_1*Tr_0*b_1*d-2*In_0*S_1*b_1+N*((a+g)*In_2+In_3), ((Tr_1*d+In_1)*b_0+b_1*(Tr_0*d+In_0))*S_0+S_1*(Tr_0*d+In_0)*b_0+N*S_2, b_2, -Tr_4+7878684013845925923520645929070239909036396792955406416139221387807839163067450930209828273643922704451797884506507398060905724506738119235258652775/20886710157231560147868050585210119183, -In_4*g+Tr_4*nu+Tr_5, ((-Tr_0*b_3-3*Tr_1*b_2-3*Tr_2*b_1-Tr_3*b_0)*S_0+(-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*b_0+(-3*Tr_0*b_2-6*Tr_1*b_1)*S_1-3*Tr_0*S_2*b_1)*d+(-In_0*b_3-3*In_1*b_2-3*In_2*b_1-In_3*b_0)*S_0+(-In_0*S_3-3*In_1*S_2-3*In_2*S_1)*b_0+(-3*In_0*b_2-6*In_1*b_1)*S_1-3*In_0*S_2*b_1+N*((a+g)*In_3+In_4), ((Tr_2*d+In_2)*b_0+(Tr_0*b_2+2*Tr_1*b_1)*d+In_0*b_2+2*In_1*b_1)*S_0+((2*S_1*Tr_1+S_2*Tr_0)*d+S_2*In_0+2*S_1*In_1)*b_0+2*S_1*Tr_0*b_1*d+2*In_0*S_1*b_1+N*S_3, b_3, -Tr_5-124506243857727683073470113596801301591591684894141699835285112911443808544901809074102550730647231351253422741988448688043130591470961319581105479421970419269520982116879769593674139108505/57521038744715258361544119804851604200383884613441, -In_5*g+Tr_5*nu+Tr_6, ((-Tr_0*b_4-4*Tr_1*b_3-6*Tr_2*b_2-4*Tr_3*b_1-Tr_4*b_0)*S_0+(-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*b_0+(-4*Tr_0*b_3-12*Tr_1*b_2-12*Tr_2*b_1)*S_1+(-12*S_2*Tr_1-4*S_3*Tr_0)*b_1-6*b_2*Tr_0*S_2)*d+(-In_0*b_4-4*In_1*b_3-6*In_2*b_2-4*In_3*b_1-In_4*b_0)*S_0+(-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1)*b_0+(-4*In_0*b_3-12*In_1*b_2-12*In_2*b_1)*S_1+(-4*In_0*S_3-12*In_1*S_2)*b_1-6*b_2*In_0*S_2+((a+g)*In_4+In_5)*N, ((Tr_0*b_3+3*Tr_1*b_2+3*Tr_2*b_1+Tr_3*b_0)*S_0+(3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*b_0+(3*Tr_0*b_2+6*Tr_1*b_1)*S_1+3*Tr_0*S_2*b_1)*d+(In_0*b_3+3*In_1*b_2+3*In_2*b_1+In_3*b_0)*S_0+(In_0*S_3+3*In_1*S_2+3*In_2*S_1)*b_0+(3*In_0*b_2+6*In_1*b_1)*S_1+3*In_0*S_2*b_1+N*S_4, b_4, -Tr_6+281080385148972428240436384088298666195348399128684395386544455996834076875384683519634664262485098546664898376585117735992463821615722882960556082361115626963265927349011318896186885552636655183974219686775635514272564045418205/22630042010746383872527102247293515467571649367959998220570201, -In_6*g+Tr_6*nu+Tr_7, ((-Tr_0*b_5-5*Tr_1*b_4-10*Tr_2*b_3-10*Tr_3*b_2-5*Tr_4*b_1-Tr_5*b_0)*S_0+(-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*b_0+(-5*Tr_0*b_4-20*Tr_1*b_3-30*Tr_2*b_2-20*Tr_3*b_1)*S_1+(-30*S_2*Tr_2-20*S_3*Tr_1-5*S_4*Tr_0)*b_1+(-10*Tr_0*b_3-30*Tr_1*b_2)*S_2-10*b_2*Tr_0*S_3)*d+(-In_0*b_5-5*In_1*b_4-10*In_2*b_3-10*In_3*b_2-5*In_4*b_1-In_5*b_0)*S_0+(-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1)*b_0+(-5*In_0*b_4-20*In_1*b_3-30*In_2*b_2-20*In_3*b_1)*S_1+(-5*In_0*S_4-20*In_1*S_3-30*In_2*S_2)*b_1+(-10*In_0*b_3-30*In_1*b_2)*S_2-10*In_0*S_3*b_2+N*((a+g)*In_5+In_6), ((Tr_0*b_4+4*Tr_1*b_3+6*Tr_2*b_2+4*Tr_3*b_1+Tr_4*b_0)*S_0+(4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*b_0+(4*Tr_0*b_3+12*Tr_1*b_2+12*Tr_2*b_1)*S_1+(12*S_2*Tr_1+4*S_3*Tr_0)*b_1+6*b_2*Tr_0*S_2)*d+(In_0*b_4+4*In_1*b_3+6*In_2*b_2+4*In_3*b_1+In_4*b_0)*S_0+(In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1)*b_0+(4*In_0*b_3+12*In_1*b_2+12*In_2*b_1)*S_1+(4*In_0*S_3+12*In_1*S_2)*b_1+6*b_2*In_0*S_2+N*S_5, b_5, -Tr_7-13325675804014921669260572317817955192541547499843631518257190167921130840758222958955188566407850063567475164697411159159815212168342859401745754594612660237106207977851944608958519249755545349384854228034431444068988714943373422403236130575425515132124579706019741265/186966283368085724585964723518518566038346933312124127222781383220456251781, -b_1, -b_2, -b_3, -b_4, -b_5, 237192631520637902028639640469025429938198712429523171375494725772017205134975705073441336611020175841889001989763396870406884061817825/20886710157231560147868050585210119183-In_3, N*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, b_5, S_5, Tr_4, In_4, b_4, S_4, Tr_3, In_3, b_3, S_3, Tr_2, In_2, b_2, S_2, Tr_1, In_1, b_1, S_1, Tr_0, In_0, b_0, S_0, z_aux, w_aux, N, a, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [In_3, b] 262
# 130, [1 = 5, -S_0 = 1, In_2*S_2 = 1, S_4*Tr_0*d = 1, -5*S_4*Tr_1*d = 1, -In_0*g = 1, -S_5*Tr_0*d = 1, z_aux*N = 1, N*In_5 = 1, -6*In_2*S_2 = 1, N*In_1 = 1, S_0 = 1, -3*In_2*S_1 = 1, In_2*S_1 = 1, -4*S_1 = 1, -Tr_1*d*S_0 = 1, Tr_2 = 1, -g*In_5 = 1, -10*S_3*Tr_2*d = 1, nu*Tr_1 = 1, N*In_2 = 1, -g*In_2 = 1, -In_2*S_0 = 1, N*S_4 = 1, -Tr_7 = 1, nu*Tr_5 = 1, -S_0*Tr_4*d = 1, -S_3*In_0 = 1, In_5*N*g = 1, S_1*Tr_0*d = 1, -2*S_1*In_1 = 1, -10*S_2 = 1, In_2*S_0 = 1, In_0*N*g = 1, S_0*Tr_4*d = 1, S_0*Tr_3*d = 1, N*In_4 = 1, -10*S_2*Tr_3*d = 1, -S_2*In_0 = 1, N*g = 1, nu*Tr_3 = 1, -S_1*Tr_0*d = 1, S_3*Tr_1*d = 1, In_0*S_0 = 1, -In_0*S_0 = 1, In_1*S_3 = 1, S_3*In_0 = 1, -Tr_0*S_0*d = 1, -S_0*In_1 = 1, -Tr_5 = 1, nu*Tr_2 = 1, nu*Tr_6 = 1, -In_5*S_0 = 1, -S_0*Tr_5*d = 1, Tr_4 = 1, -In_4*S_0 = 1, Tr_0*S_0*d = 1, In_1*N*a = 1, N*S_3 = 1, Tr_6 = 1, -3*In_1*S_2 = 1, -S_1*In_0 = 1, -5*S_1*Tr_4*d = 1, nu*Tr_4 = 1, -S_2*Tr_0*d = 1, -4*S_3*Tr_1*d = 1, Tr_3 = 1, Tr_1 = 1, -3*S_1*Tr_2*d = 1, S_2*Tr_2*d = 1, -g*In_6 = 1, In_0*S_4 = 1, In_4*N*g = 1, S_0*Tr_2*d = 1, N*a = 1, Tr_1*d*S_0 = 1, -2*S_1*Tr_1*d = 1, In_2*N*a = 1, -5*In_4*S_1 = 1, In_1*N*g = 1, -Tr_2 = 1, -3*S_2*Tr_1*d = 1, -6*S_2*Tr_2*d = 1, In_4*N*a = 1, -S_3*Tr_0*d = 1, S_1*Tr_2*d = 1, -Tr_1 = 1, -S_0*Tr_2*d = 1, -Tr_4 = 1, Tr_0*nu = 1, In_5*N*a = 1, -4*S_1*Tr_3*d = 1, S_1*Tr_1*d = 1, N = 1, N*S_2 = 1, S_1*In_1 = 1, S_1 = 1, -10*In_2*S_3 = 1, -Tr_3 = 1, -4*In_1*S_3 = 1, In_4*S_0 = 1, -g*In_4 = 1, -g = 1, -2638160105128370917320815258204546212760429671424506002980329533587986305429272675295666893652158029195494/82303096410101184534590041 = 1, In_0*N*a = 1, -Tr_0 = 1, S_3*Tr_0*d = 1, -Tr_6 = 1, S_1*In_0 = 1, -5*In_1*S_4 = 1, -11789944372886487387107349367984574223504359631571202251428134167414008940707182513136978547404552574493972851075718769518558642821926401190547010305669148073602172455220303556136600810458521423570255056859140853099201035935442825382265107639538555267631480548205074/1672514064053907842387051727427847486825017505927337569102457886044552712316763 = 1, -In_0*S_4 = 1, In_2*N*g = 1, Tr_7 = 1, -In_0*S_5 = 1, -g*In_1 = 1, S_0*In_1 = 1, -101822985103613257601787337477588801128388264758326686794961002882236082107892833540580618234697336904653124219695418089796838618562753784891502689466878109169917341661140182023753659206/6773799678690410489111140486807745823203599936381681 = 1, Tr_5 = 1, S_2*In_0 = 1, N*In_6 = 1, -S_0*Tr_3*d = 1, S_2*Tr_1*d = 1, In_1*S_2 = 1, -S_4*Tr_0*d = 1, S_2*Tr_0*d = 1, N*S_5 = 1, N*S_1 = 1, S_1*Tr_3*d = 1, -1 = 1]
# 134, -4.837786147
# 1
# [In_3 = 9, b = 67]
# [In_3 = [2, 2, 3, 3, 3, 3, 3, 3, 3], b = [4, 3, 4, 4, 3, 3, 4, 3, 4, 4, 4, 3, 3, 3, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3]]