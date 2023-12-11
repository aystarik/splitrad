#pragma once

#ifdef NEED_TEST_DATA

#include <stdint.h>

static uint16_t a[] = {
    6676, 15256, 25569, 32461, 33604, 35810, 42034, 46825, 56732, 56298, 49104, 44762, 44021, 43685, 45729, 47986, 45337, 40520, 27877, 20547, 17100, 15275, 12852, 14203, 12374, 9384, 5263, 6379, 11642, 10666, 5963, 713,
    10944, 13556, 21658, 29659, 29491, 33267, 41084, 44591, 51212, 51844, 47210, 43159, 40890, 37875, 41730, 45404, 41574, 32219, 25627, 22984, 19160, 17885, 14897, 13466, 13289, 11490, 11179, 10567, 7043, 5530, 3366, 6837,
    19493, 23239, 25094, 28292, 27147, 31201, 37253, 43449, 49010, 51157, 50977, 45499, 35512, 33211, 37540, 39584, 36515, 29324, 24860, 21146, 17977, 21418, 18866, 17744, 20525, 20921, 15307, 14899, 12670, 10929, 9170, 9334,
    24087, 26130, 28839, 33654, 30709, 29478, 31936, 34980, 43085, 49769, 50709, 42884, 34106, 33115, 34997, 36179, 29952, 20114, 17825, 15997, 16143, 19375, 17780, 21029, 23992, 26294, 21867, 20290, 20743, 18476, 17386, 21203,
    29731, 24644, 23908, 25756, 28221, 29128, 26793, 29521, 40310, 44839, 46538, 41421, 33981, 32231, 26676, 30831, 30347, 18744, 11026, 11848, 14587, 16976, 19089, 19122, 23388, 26320, 27395, 28889, 30162, 31797, 28033, 32616,
    32988, 22302, 19606, 19098, 24711, 22118, 21212, 24625, 36706, 41023, 40813, 35050, 32027, 28901, 26235, 28607, 29773, 15052, 6505, 9700, 14691, 16202, 22077, 21657, 24226, 27855, 32091, 36347, 38924, 42636, 42691, 36451,
    39163, 26316, 19768, 17013, 23980, 23150, 21258, 23531, 32544, 39194, 37655, 32481, 30398, 28008, 20867, 22542, 20661, 8962, 7969, 8770, 7659, 12090, 15493, 20664, 26245, 27505, 30531, 36764, 43089, 48842, 51388, 44159,
    43133, 31390, 24791, 19341, 18685, 20887, 20148, 21695, 21692, 24710, 28320, 27929, 28708, 26516, 22260, 19994, 11105, 7653, 9415, 7320, 1746, 6285, 8654, 19769, 24132, 24729, 31391, 37420, 42580, 47974, 49929, 47374,
    46547, 38022, 28141, 22140, 19142, 20225, 15096, 13349, 10313, 4394, 14390, 21562, 27352, 28208, 21285, 14545, 9069, 10727, 7771, 9489, 4355, 3049, 10748, 20473, 25675, 32940, 37735, 39905, 43342, 47745, 48142, 47436,
    46136, 39031, 30442, 31721, 25949, 18874, 19025, 14030, 10255, 3137, 2609, 15986, 22818, 22087, 16149, 10627, 9038, 9788, 9029, 10462, 10374, 12581, 17395, 24020, 24763, 34062, 42531, 41523, 41316, 44792, 46512, 48928,
    46818, 37724, 34677, 33924, 25050, 18124, 20741, 21937, 12977, 5684, 3268, 10395, 11551, 10963, 10911, 6051, 6088, 12242, 15914, 13710, 16753, 18361, 19138, 23018, 26300, 34167, 40625, 41597, 41994, 46755, 46524, 47386,
    41506, 39245, 39197, 37285, 26958, 20579, 24564, 23908, 18662, 9785, 7209, 8199, 9454, 7187, 5402, 6228, 11370, 19625, 20933, 16170, 16263, 19045, 26041, 26310, 28562, 33666, 38426, 37124, 39176, 43524, 44034, 42688,
    38444, 37489, 38445, 38775, 31268, 25316, 21646, 24152, 18194, 11834, 7580, 10903, 14385, 14909, 9824, 7865, 12279, 14480, 19304, 14870, 16110, 19151, 25453, 29217, 33546, 33874, 39249, 39149, 40980, 43171, 45105, 43081,
    41390, 36319, 36687, 37177, 32209, 29212, 26193, 26891, 20202, 15608, 13645, 14200, 15302, 12946, 11845, 11983, 7488, 11431, 21217, 22070, 20603, 24051, 23186, 29443, 32374, 37160, 43098, 44466, 43889, 44474, 49056, 46904,
    44025, 34293, 36016, 32376, 31360, 28157, 25479, 24164, 21813, 25710, 22462, 17701, 15916, 6521, 5309, 9130, 10434, 12192, 23228, 28022, 24931, 26336, 27663, 33805, 36903, 37358, 42187, 42943, 48408, 48989, 56209, 56675,
    50474, 35176, 32318, 27556, 26775, 27545, 24906, 16598, 14689, 24746, 28959, 21285, 11950, 3447, 3363, 11495, 13339, 14167, 16227, 23267, 27640, 26124, 28693, 30004, 34656, 38622, 40963, 41220, 49153, 56655, 64739, 63987,
    52686, 41472, 34404, 27073, 22493, 20828, 17577, 10806, 9614, 21619, 26196, 19920, 11792, 6964, 6136, 13118, 15131, 11900, 12225, 17026, 21804, 22731, 23160, 26023, 31084, 36540, 32298, 38927, 48582, 57558, 64447, 60091,
    49221, 38842, 32893, 28412, 19952, 13662, 10758, 11152, 12549, 17755, 19735, 13178, 6473, 4550, 8623, 9101, 13059, 10807, 7012, 12002, 17355, 19131, 24518, 29863, 30850, 36433, 37063, 40388, 47456, 53501, 53544, 56719,
    45997, 33152, 25481, 20938, 14576, 5539, 3628, 7945, 12439, 13596, 12333, 6213, 1105, 12269, 12472, 9680, 12689, 15304, 11231, 15240, 17464, 22420, 30765, 35999, 38693, 38248, 43085, 45517, 48693, 51814, 47356, 50116,
    40911, 31051, 18867, 15760, 7790, 3065, 2874, 6169, 16534, 17062, 8343, 7631, 9489, 18754, 23987, 16748, 18239, 14012, 15796, 18658, 24128, 33072, 33527, 36436, 38567, 37524, 38629, 41919, 45094, 46652, 47178, 44874,
    33379, 23197, 10536, 5179, 7059, 7922, 8406, 10844, 21309, 25611, 21890, 21990, 25248, 33261, 36849, 30207, 23256, 20409, 24606, 25249, 31161, 35627, 32930, 32418, 33381, 35551, 33259, 35874, 42235, 43587, 44996, 41067,
    24399, 14582, 5257, 8510, 11538, 16135, 15793, 23392, 28311, 26651, 27797, 29940, 34722, 38265, 37570, 36317, 31543, 26484, 27538, 28136, 30840, 32888, 31985, 29516, 26094, 28044, 24670, 27023, 36755, 42947, 44106, 37526,
    19249, 8744, 8500, 17960, 16569, 21905, 24379, 29817, 34275, 34526, 36365, 38292, 39341, 36859, 36056, 40397, 34142, 23153, 23488, 26364, 32662, 31434, 28470, 27360, 23139, 20617, 20544, 22290, 24703, 31140, 36561, 28415,
    17897, 9266, 18721, 28903, 30595, 33101, 40365, 41100, 44968, 45151, 45338, 47719, 49176, 41287, 34977, 33422, 27241, 14291, 16179, 22492, 28560, 31529, 27221, 26664, 25304, 19966, 18501, 17095, 18303, 21860, 24584, 22598,
    20178, 23074, 29721, 37353, 39091, 41038, 49452, 51877, 51653, 51674, 50647, 52093, 53628, 43606, 34412, 28095, 15586, 2293, 9986, 19723, 26977, 32194, 26353, 24772, 22605, 18898, 17171, 12473, 10283, 13608, 17535, 16929,
    23848, 32789, 41827, 45470, 45564, 45126, 52017, 54522, 51682, 49995, 50043, 52313, 53200, 40139, 31859, 22883, 13512, 3509, 7099, 15453, 20935, 21592, 20199, 19295, 15008, 15430, 14455, 13954, 8608, 5545, 12598, 15495,
    20255, 31757, 43065, 45089, 44746, 46648, 52495, 48586, 46256, 47407, 52621, 53221, 49563, 41533, 37572, 27529, 14690, 8550, 3876, 14469, 17588, 13644, 14506, 12315, 6009, 7639, 11667, 14894, 6852, 1461, 6264, 10629,
    18579, 29788, 42245, 47352, 48098, 46928, 48200, 43945, 44886, 51393, 55705, 59589, 56500, 52756, 45896, 31850, 21247, 12433, 1065, 6724, 10140, 11441, 10468, 5291, 2751, 4624, 7732, 5949, 3739, 0, 5845, 12921,
    22179, 29766, 39949, 48106, 51121, 51931, 46649, 46436, 47162, 53042, 58884, 65535, 61713, 58164, 48569, 38707, 31959, 18308, 2719, 2145, 3274, 8957, 11275, 8925, 6945, 6083, 6045, 8068, 10341, 2397, 4335, 16428,
    21859, 28890, 33177, 39268, 49162, 50280, 47134, 49953, 46675, 47987, 54022, 59454, 60714, 56590, 54167, 40206, 32015, 23595, 7192, 6039, 7727, 15781, 14706, 18649, 18153, 14350, 11558, 10965, 8619, 684, 7777, 16870,
    19404, 28007, 32128, 34378, 46429, 48556, 49871, 50698, 50155, 52950, 48017, 48349, 53239, 53877, 54837, 45756, 35718, 27087, 16785, 9645, 11786, 14331, 16032, 20498, 17177, 13721, 9772, 10516, 10620, 4429, 4603, 10210,
    14760, 25181, 30962, 34048, 41506, 44581, 48497, 51465, 57616, 57314, 52385, 44419, 45959, 45177, 45728, 47928, 43452, 36111, 26528, 16748, 13512, 12968, 13611, 17604, 12516, 7083, 7255, 14964, 15957, 9777, 5754, 6269,
};
static uint16_t b[] = {
    41636, 44907, 43477, 45893, 41648, 39514, 35671, 38693, 37101, 30974, 26704, 25158, 23082, 15911, 11891, 8749, 15172, 15979, 17216, 12718, 10705, 10912, 15793, 19370, 16198, 17713, 22363, 24994, 30646, 34711, 35900, 41932,
    44315, 44842, 47243, 50656, 47196, 39275, 34870, 36776, 33980, 29510, 25954, 26401, 25000, 19855, 18439, 16146, 15470, 15238, 9849, 10475, 10928, 5812, 13142, 24689, 24319, 24371, 25739, 25563, 33615, 33231, 38711, 42774,
    42730, 51888, 53810, 61302, 58177, 41584, 34447, 33136, 30369, 29950, 24994, 23395, 20320, 21253, 28908, 23525, 18908, 13089, 3787, 3475, 9236, 11287, 14983, 21952, 26890, 25085, 24926, 29885, 36796, 37297, 39602, 41320,
    42944, 51244, 60886, 65535, 60797, 49729, 35372, 33367, 25383, 24901, 26234, 22274, 12514, 13596, 26157, 29218, 19992, 11088, 4047, 8181, 14856, 13684, 13370, 16224, 23989, 27221, 26930, 28189, 28170, 34935, 39078, 37840,
    41302, 48911, 57385, 63264, 58027, 48658, 39190, 34364, 26216, 21311, 19633, 13175, 10407, 14058, 21436, 22229, 17035, 10895, 6109, 7774, 13824, 14807, 9375, 12834, 16140, 23425, 21253, 23256, 26844, 32377, 34847, 32576,
    42811, 51525, 52018, 50558, 55050, 45966, 35101, 28694, 24459, 16278, 10850, 9168, 11647, 12858, 18239, 17065, 8801, 2284, 5937, 10293, 10112, 14631, 11266, 6952, 13483, 18556, 19798, 28624, 32765, 33178, 38055, 40278,
    44842, 47934, 49030, 48325, 49369, 42179, 27911, 21154, 18542, 11115, 2510, 875, 8752, 14291, 11437, 9098, 4650, 4638, 17310, 12308, 11777, 16287, 15022, 14432, 18542, 20140, 27319, 32711, 38669, 41129, 39840, 43229,
    43760, 43767, 46691, 46363, 41524, 38935, 24494, 15596, 12193, 5456, 3992, 4763, 9491, 21988, 20009, 12465, 13361, 15100, 26702, 27191, 21142, 19067, 14065, 18036, 21467, 30081, 35795, 34219, 36854, 35310, 36841, 36489,
    35786, 42518, 44746, 45005, 38525, 27750, 18047, 6290, 5925, 10132, 10862, 12094, 17915, 25150, 26998, 24326, 26022, 32747, 38331, 38354, 30198, 25105, 24991, 27653, 28270, 32301, 35345, 32511, 31479, 31729, 32878, 30872,
    27409, 34655, 42328, 42000, 32823, 19450, 9360, 7110, 13103, 13303, 20103, 18827, 28492, 28813, 28363, 30238, 34916, 38653, 38477, 37329, 37021, 31018, 26612, 27104, 28552, 33688, 32995, 31096, 27967, 22743, 23844, 22594,
    21852, 22483, 30489, 32137, 23682, 15160, 6769, 14945, 22767, 21670, 26294, 30009, 34456, 37411, 38448, 40624, 40982, 41479, 38069, 37593, 37606, 29488, 19454, 22237, 27617, 32455, 30368, 27075, 25068, 23833, 19061, 20777,
    15305, 18604, 21229, 20580, 21394, 17976, 13144, 25275, 33217, 36168, 38054, 45069, 46363, 47611, 48537, 46878, 51073, 50221, 40292, 35698, 31020, 19233, 9723, 15679, 22755, 29240, 33046, 26362, 24988, 24150, 18283, 17432,
    12328, 9912, 11787, 17886, 18330, 24066, 28740, 35566, 41800, 42297, 44935, 52198, 52710, 51999, 51437, 50771, 53660, 52244, 41069, 32073, 24588, 10323, 1778, 11753, 18993, 28260, 28363, 23807, 23881, 20274, 18675, 15637,
    14786, 6771, 5835, 10438, 14920, 26950, 35355, 44909, 46492, 44512, 47722, 54340, 52822, 49521, 48226, 49386, 54480, 49450, 36428, 31734, 20387, 12553, 4336, 8163, 15724, 21886, 17469, 19244, 16751, 10999, 13771, 13783,
    10880, 4239, 572, 7971, 11880, 21879, 33468, 43647, 43369, 46095, 47031, 51607, 45645, 45015, 48777, 54835, 54335, 49456, 42533, 40100, 22633, 14463, 7239, 3647, 15066, 16911, 12833, 11286, 8632, 5693, 7153, 11475,
    5261, 5685, 2008, 7575, 16401, 22850, 33583, 43653, 48390, 50221, 48587, 45942, 43197, 47610, 53251, 59569, 62253, 58327, 54661, 43108, 31938, 23766, 11484, 443, 5176, 7479, 12342, 10555, 4996, 4276, 7540, 5076,
    13077, 11095, 0, 8166, 19710, 24892, 29470, 39925, 47171, 51450, 49155, 46437, 47460, 47058, 52327, 60927, 63527, 61235, 56161, 47220, 37459, 29839, 16965, 1794, 3016, 5601, 11324, 12440, 12321, 9793, 8182, 10606,
    8446, 7532, 3095, 10608, 16781, 22001, 31006, 31922, 39435, 50513, 47106, 49123, 50383, 47727, 48745, 51740, 58671, 57520, 57411, 53181, 39173, 31661, 21374, 8840, 7009, 10858, 18890, 16046, 21860, 19161, 17715, 12338,
    14529, 13471, 4608, 7295, 11945, 22268, 27939, 32086, 35348, 46500, 49319, 50052, 50556, 53125, 54359, 46075, 46435, 50597, 52938, 50826, 46744, 36110, 27066, 19415, 11523, 12512, 13302, 17537, 19285, 12813, 11619, 12142,
    15985, 15852, 10401, 5919, 3526, 13625, 23488, 30190, 34626, 39480, 43888, 47654, 52345, 59971, 56251, 53032, 44228, 46298, 45322, 45117, 47719, 45368, 37129, 24700, 15992, 15595, 15256, 14189, 17423, 10434, 8146, 8184,
    7864, 10505, 9155, 3273, 2134, 7632, 14569, 27819, 34510, 31124, 35655, 42332, 46042, 56649, 53494, 45147, 43257, 42596, 42584, 46518, 47048, 44158, 37295, 25314, 21311, 18664, 15538, 12485, 13984, 12040, 9483, 4770,
    14038, 6849, 6870, 7020, 11600, 14783, 18897, 24447, 26802, 28656, 35446, 43071, 45931, 49762, 49995, 46291, 43704, 38386, 37287, 39971, 42312, 37528, 30045, 26135, 21799, 19627, 19515, 15722, 16035, 14843, 13107, 13396,
    17668, 14510, 12212, 10740, 12280, 24580, 24330, 28121, 30181, 27300, 31141, 36557, 41666, 47729, 52757, 50913, 43056, 31903, 33026, 38562, 37803, 33050, 23695, 24313, 18813, 18168, 21859, 18819, 20463, 22920, 22304, 17389,
    23590, 24053, 20716, 22484, 28605, 24626, 26054, 30328, 32446, 30080, 31720, 30162, 35296, 42378, 50366, 48002, 38530, 33492, 31618, 32566, 32546, 26646, 17118, 15128, 14411, 19383, 18575, 18142, 22278, 26463, 26151, 22634,
    30948, 34257, 34958, 31715, 33403, 29093, 23830, 21236, 22993, 26456, 26135, 26548, 32582, 41360, 42248, 44871, 37148, 33849, 28878, 26345, 31620, 30535, 14979, 10237, 12340, 15492, 17825, 21708, 19187, 23287, 27088, 29892,
    37528, 43332, 46926, 45112, 38679, 32612, 22826, 19246, 20256, 24821, 20198, 20123, 25575, 36523, 40044, 37972, 33245, 30256, 27864, 23285, 28474, 24533, 11684, 5756, 13715, 14366, 15874, 22155, 23829, 28272, 29483, 32914,
    38664, 43432, 51084, 50987, 46163, 37858, 24524, 20355, 18136, 24023, 22363, 21529, 24414, 32290, 36399, 33829, 29093, 29002, 24440, 18558, 21957, 14104, 9353, 11151, 6562, 6042, 11446, 13890, 22214, 24259, 28146, 30971,
    39405, 42810, 47920, 48921, 48697, 42976, 32439, 26385, 17781, 19328, 19894, 17080, 18601, 16993, 19345, 23921, 23878, 28743, 27600, 22675, 15335, 8945, 9289, 8188, 7349, 272, 6709, 10369, 21750, 24773, 28541, 36399,
    41436, 45389, 48412, 45965, 48511, 45950, 37209, 27825, 25839, 20814, 19536, 14491, 11875, 6964, 2967, 12927, 21116, 27592, 26081, 18138, 12392, 10032, 11166, 10764, 9583, 6516, 8190, 16052, 22197, 28515, 37003, 40772,
    44209, 42656, 44924, 47597, 50090, 45269, 34616, 31053, 32745, 24145, 18025, 20189, 15139, 9808, 3543, 2573, 15454, 20352, 18903, 13911, 9697, 9400, 9594, 12348, 10642, 12718, 16170, 19771, 23542, 26787, 35333, 41377,
    40420, 42359, 47017, 46617, 47023, 44428, 36118, 36440, 32863, 22546, 19458, 22924, 21971, 12684, 5030, 6346, 9761, 10197, 6585, 8766, 6453, 9835, 16721, 19534, 14195, 20324, 21487, 22026, 23721, 29087, 34218, 39793,
    36892, 40007, 43872, 44477, 40424, 40380, 37416, 39196, 37530, 27024, 21286, 24140, 23462, 17948, 10632, 8773, 9743, 11303, 8769, 6292, 9147, 14834, 21222, 21167, 15205, 17224, 21847, 27908, 28845, 34174, 34932, 39762,
};

#endif
