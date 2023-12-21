import os
import tempfile

from socialgene.base.socialgene import SocialGene

from .test_mibig_gbk_parser import PROTEIN_DICT

DIRECTORY_OF_THIS_FILE = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(DIRECTORY_OF_THIS_FILE)

FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "test_genomes")


def test_fasta_file_parse():
    # just create the fasta from the genbank file
    with tempfile.NamedTemporaryFile(suffix=".faa") as fp:
        sg_object = SocialGene()
        gbk_path = os.path.join(FIXTURE_DIR, "lagriamide_mibig_bgc0001946.gbk")
        sg_object.parse(gbk_path)
        sg_object.write_fasta(outpath=fp.name)
        fasta_object = SocialGene()
        fasta_object.parse(fp.name)

    protein_parse_results = {
        k: [v.description, v.external_id, v.domains, v.sequence]
        for k, v in fasta_object.proteins.items()
    }
    assert protein_parse_results == {
        "Tdc2m3PRLsyEzjwyux6BF4arDy2mQ_Bl": [
            "AXA20086.1",
            "AXA20086.1",
            set(),
            "MQEYCRLRRKLTLELSPEDAADVAQEAFERTLRYMRKHDGRVASPVGLLVRIALNLQIDRGRRRKHLPTALDESWEHPRWDITPEDEVVGRQSVTQLVETLDKLAPRRREAFVLCRLHGLTYQDAAKKMGIRPSVVREYLVDAVRACRDSVDWAVSRVKCNTPLVNGLELSRSG",
        ],
        "ptq1NGhBcUp3TIEqvAUxnnp4LOKwINvn": [
            "AXA20087.1",
            "AXA20087.1",
            set(),
            "MAPTGIILGDTNGLRISLVCNPFVSPKIMPVGAIASSRSIYATIPLAMDGPVYWILWTDNPSNR",
        ],
        "-l7xLyFZbiZENPLq_GML8JyTRF1Srawr": [
            "AXA20088.1",
            "AXA20088.1",
            set(),
            "MSEYIHKSHNVTVLTYHVVLWRGIEKQSLMTGLAKY",
        ],
        "T_DzOorDp3ROhRRBtuXP3xyAPorpTVD0": [
            "AXA20089.1",
            "AXA20089.1",
            set(),
            "MSDLPPRFAKSAAVERRFSQTKQGCNAGFARAMPGLVSFFTPQFADMNKLIDYWRATVACGRMGEGMSERGRHTEFSVVLLRNRDSAPVFGVVLRA",
        ],
        "AStsOnOU5ZWxURs9PrTiWjddkuQXfanl": [
            "AXA20090.1",
            "AXA20090.1",
            set(),
            "MQSVNPIGETKKALRLLQELVGRELSLPADEVPVDTDYLNMGFSSRSLIGLIQQLSSLLRTKLNPSVLFDYRTLTEFADYLAEHHAQWLTAIEREQPAERNQEAGDETSVAVPLSEAQRGLWFLQKNHPWMAAYNVPLCLRLSPRVQRERMRQACAWLPRRWPVLGASVQRADGRLVMQTQPARKLTWQEHHAEDWSEAQRLAWLGDRLAEPFDVDNGPLLRAYWLGGEPDGASRLLLVIHHLVIDAVSVGVLLAGLRKTYADLEGWRDLSGAVDDSAYGAFVAAEAERLAGAEGFARLAYWREQLADVPGSLGLPLDRARGVTPSFKGCTLRRELQAALGESLRAYTERHRVYPSTLLLAVFQGLLSRHAGRDDVVVGMAMDERDAASAGLVGMFVNMLPMRARGLGRRGFVEDMQALQRQLVDAMAQAYPFPALVRELGLSGSDASPLFQAAFLYQDTLDIDVLNGVDDWVWEEALYQEGEYELVLEVRRRLKGYALYFKYDPTLWDKSTIERWLAHYLRLLEGVLAAPQKRLGEHELRGEHERARLAALQGEVRDWALTPVSILFERQAIANAQAWAVSDDQQRWCYAELAAHSAAIAQRLHVQGIGTGSIVGVCQGRSPWLLASLLGIWQAGAAYIPLDPAYPVERLRYMLEDSGASAVLSDTSHLVLVQALAGMLPVWAADAAAPSTSAPSFPPPQPEDLAYVLYTSGSTGRPKGVRISHGALSNFLQSMAEKPGLTAGDRLLAITTISFDIAGLELYLPLIGGGECVLCPEEIVRDARRLKAEVEHVRPTLMQATPATWSMLFHAGWRNAEELKVLCGGEALPARLKQRFDEIGTTVWNLYGPTETTIWSTLAKLDAEDTSIHIGQPIANTQAYVLDQEGREQPIGIEGELYLGGAGLAQGYHGQPERTAQAFIDHPLGRLYKTGDLARWRADGQLEYRGRSDQQVKVNGHRVEPGEIEAVLEQSGLVKQAAIVLREGAHGSQLAAWCVPTKVTQGDTWLDPVQIQVLQAQLRDRVPAYMQPSIWLGTAALPKTLNGKIDRQVLSARALPEQREEKAPPSAVTRSAARRALESRLQALWAQVLERSTVSRDERFMEAGGNSVAAVLLAERIQAEFGRAFGVAQVFAYTSVAAQAAYLDASDLTFEALANEPIATQVIVAEPERDASGVAEDALAIIGIACHFPGAEDHRAFWNNLRAGHDSGKLFSPEALRAAGVPERLIADPHYVPIRYGIEGKAEFDADFFNLSPRAASLMDPQYRLLLQQAWAAIEDAGYTPEVIPDTGVFMSASFSAYQARLQDPSAVEAGDRYVAWLMSQGGSLPTLISYHLGLTGPSLFIQTNCSSFLSALAAARSSLLARESRLALVGAATLFTEDSKGYLHEPGLNFSSDGHCKTFDASADGMVSGEGCGVILLKQAREAVADGDHIYALLRGIAVNNDGADKAGFYAPSVRGQSEVIEQALRQARIDPGQIGYVEAHGTGTRLGDPIEVTALSETYRRHTQATQYCAIGSVKSNLGHLDTAAGLAGCIKLALSLQHGEIPPTLHYQQPNPAIDFAASPFYVAEHLQAWPPGPRLAGLSAFGIGGTNTHAILEAYLEESVAARVDGAQLIVLSARTQERLHAVERQLLDHLDSSASLPSLRDLAYTLQVGRRGMTHRLAFVVEDVSALRRQLSDHLAGRKVGHEGDCDRGHAVLHAFENDGDSARLFAQWAAAGKLDKLAALWVQGLSLDWALLHGGQRPRRVSLPTYPFERTRHWLEVPLAASAVSPESPELATERSWDGMSYLPRWVATPATAPDVEEVPAPASLLIVAPEASARADELFEWCTRRWPSAALRLIRLGRENLWLGEAERVCDSFDDGQALTEALREGTWPDSVVFLAEEAGSTLNWPGHPESAEAQLLRVQQALSQSQPAQRIEFYLVTIEANASDALTGGGLSGLAYALAQGDHRLRLRHLSLDADVWSASSWWQLWAEPASDRGERVRFSGGERWRQRFDRLNWGSLRDGGLKQGGVYLIAGGSGTLGIAISRYLIERYRARVIWLGRSRADAPELVARRQVLEAGALLPGYVQADLTDATAVHEAVARARQIHGTIDGAIFSAMYRGADAPAERWSPTALAGAVAVKALGARHFYQALLGESLDFLCYYSSVQSFAFLSARNSAAYAVGVAAADRYVQLIRAQSPFRVGLIHWGYWQDTLAGTALEQDLARHFTGIHADEACVFFERFVAALGQGMLDQAICLKASDTVRDVMPSTVDDVVVLTTAGTPSFLDGFDAAALRSTVMPPDWTELDRRLRGLLCAQLRYLGLFDQAGVSWESEQLRRKLGVIDDYRRWWDECCTEMLEKEGWIRTRKGRVELLQALSLEQADALWQDWEHDKADYLKDPQLRAALLLVEDCLRHLPAVLLGRTSATSVIFADGSMSKVAGLYQGTAWTDSFNHQVADTVEVYVRHRLAADPLARLRVLEVGAGTGGTTAKVLPRLAAIGAPMLEYCYTDLSEAFLAQARERHLANYPYLRTCRCDIERPLAEQGIVPGSYDIVIATNCLHATHDIRATLRHVQAALRCHGVLIANEGVSKSLLGTLTFGLLEGWWLYDDPQLRIPGSPLLDSAHWRALLDEAGMRPVCLDGPGRELQQVWVAQSRGLIRPGGFASASATRLPAPSSVPAKAGVAKPVVSASPVSSSEVVPYAAIAAEIRACLAETLKREAAGLADETAFSDYGIDSILSVALVKRLNQRLGVQLGQAVIYDYSTIAALSRHVLERRDGQAAMVPRQFGVVAASESTVVEVPVLAPVTVSEAVPPAAEAARPTPIAVIGMALQVPDAEDADTFWANLLSGHDSIRELPEAYGRPAHGASPRGGALQGRDYFDAAFFGLSNEEAAGMSPAQRLVLEEGWKALEDAGYDPRSLRGSRTAVFVGAEPSGFFQGSFSGSSEAIIASRLSYLLDLKGPALVVNTGCSSGAMAIHLACESLRRGEIQLALSGGVASALSVEGLRHLADAGMLSPQGHCLSFEATGDGMVLSEAVGMLVLKRLDQAIADGDAIHGVIAASGSNQDGTSNGITAPNGRAQEDLLGEIWARHGIAPEHISHFEVHGTGTSLGDAVEGNAISRAFGRVTARRGFCVIGTSKTHIGHTGAASGVVGLIKLLLSLRHRQLPGLLHFERANPLIDWNASALRLPEATTAWEGEPDLPRYAGLNAFGHSGTNVHMVVREYGPGEGDQRGPQLLECAQEVLLPLSAASALQLARSARRLLDFLDEEAGHGQSRPSLAELAYTCQTGRSALVERAVLKAGDREQLRALLLALAEQRPMAGLWRGSVDPETTLVAAANRDGLDELAESWIRGAEVDWQKLYGPVRPRRCHLPAYPFDRRHFPWHLAASVPSPVSRHRVPVATEPIATLTTVTAMPAWRFVLAATAEAGAEPRAQATQWLCRWLAMRLQRPVQALNPQLSYRELGLTSLGLVALSEELSRLLGVVVLPSLLFEYPSIASLAAHLAEQHAALLSRVQSIPLTGADSPSVGTPAEASVPDRVLSVLQAYRNGALNHAQARTLLAETTP",
        ],
        "IsCrCflKZgA6ghoHxXclbsOix0bbDkwZ": [
            "AXA20091.1",
            "AXA20091.1",
            set(),
            "MNAFMTLDQMIALYKAGKVSVDQVKEAFARLPDATADTGIVTDVELSETQRGLWSLQKAYPWLAAYNVPLCLRLPADLDRERFLRACAGLLERWPVLGASVEQQDGPLRLRMASVSGLSLEQDDGLAWSESERLDWLRERIEQPFDLAKGPLLRAHWLEGSGGGESLFLLVVHHLVVDGASVGLLLAGLHEAYRALGNGEAMPSASGSDGYLGFVQAERVRLQGDQAARRLAYWREQMADAPSSLGLPLDRPRSATPSFKGRTLRRELLSALGDSLNAYTEYHGVYPSTLLLAVFQGLLSRHAGCDDVVVGMAIDERDAASAGLVGMFVNMLPMRARGLGQRGFGEDVQALQRQLVDAMAQAYPFPALVRELGLSGGEASPLFQAAFLYQDMHDMLDTEALAEVSKWTWEEALYQEGEYELVLEIRRRAQGYVLYFKYDPMLWDESTIERWLAHYLHLLEGVLADPRKRLGEHELRGEHERAWLAAWQGEVRDWALTPVPALFEHQVAANAQAWAVSDDQQRWCYAELAAYGTAIAQRLHAQGIGSGSIVGVCQGRSPWLLASLLGIWRAGAAYVPLDPAYPVERLRYMLEDSGASAVLSDTLHLAQVQALAGILPVWAADAIGPLTSVPLPTLQPEDLAYVLYTSGSTGRPKGVRISHGALSNFLQSMAEQPGLMAGDRLLAITTISFDIAGLELYLPLIGGGECVLCPAEIARDGRRLKAEVERVRPTLLQATPATWSMLFHAGWHNAERLKVLCGGEALPVRLKQRFDEIGTTVWNLYGPTETTIWSTLAKLDAEDTSIHIGQPIANTQAYVLDQEGREQPIGIEGELYLGGAGLAQGYHGQPERTAQAFIDHSLGRLYKTGDLARWRADGQLEHRGRSDQQVKVNGHRVEPGEIEAVLEQSGLVKQAAIVLREGAHGSQLAAWCVPTKGTHEDTWLDPVQIQVLQAQLRDRVPAYMQPSIWLGTAVLPQTPNGKVDRRVLCARALPAQREANTQPPSATPRNAVRQELESRLQALWIEVLERPMVDRDERFMEAGGNSVAAVLLAERIQAEFGRVFGVAQVFAHTSVAAQAAYLDEVQPVSTPGETTDSSAPKPDTPETVLAEDTLAIIGIACHFPGAEDHRAFWENLRAGHDSGKLFSPEALREAGVPEPLIADPQYVPVHYGIEGKAEFDAEFFNLPPRTVTRMDPPFRLLLQHAWAAIEDAGYTPEAIPDTAVYMATGGPLRLAELTEPGSPGDSDDYVNWLLAQPGTVATMVSYQLGLRGPSYAVHANCSSSLVGMHAAAQSLRSGEVRTALVGAASLFGDDSLGYRYEPGLNFSSDGHCKPFDLRADGLVAGEGVAVVLLKRAREAVADGDHIYALLRSVAVNNDGADKAGFYAPSVRGQSEVIEQALRQACIDPDQIGYVEAHGTGTRLGDPIEVTALSETYRRHTQATQYCAIGSVKSNLGHLDTAAGLAGCIKLALSLEHGEIPPTLHYQQPNPAIDFAASPFYVAEHLQAWPPGPRLAGLSAFGIGGTNTHAILEAWPPATARLLQSEPAPGVAQVFPLSAMQADRLPVYARRLATFLRGPYAASLRLADIAYTLQTGRRSMRSRCAFVAETLDQLLAALDDCADTVDTQDAAAASPPAPADGWQLQRDAHGLAAAWQQGEPIDWNALQAVAPWSARRVSLPTYPFSPKRWPRTAARVAASTQAAVLHPLLHRNTSDFAAYRYSARFGGDEFFLADHVMGGHKVLPGVAHLEMACAAFAQAIAAPQAALELRDVVWIRPVGVDAPLDVHTVLRPQADGSADFEICSQPDAAGERIVHSQGRVMPVDTVESEAELLALDTLRQQLAQAHRDAASYYRDFEQGGASYGPAFRALEQLWLGDGQALARLVLPAARHEGAERYGLHPSLLDAAVQAAFIGINALRAAAGVASSEQGGSLPFELKRVQLLVPCEPVMWAWVRYTQGIGAGERVQRMDVDLCDEQGRVRVRLHGLARLAALVPKLAAPSLQLPQWEDAPLVYGADAPPAYAQRLLVLCGLPAHGLEAESGADLVERLEAGTGLDSVQDYALRLLARVQVFMRDKPRGRCLLQVVIPAKGEMVSLAALAGLVRCLRLEHPGVSAQLIAVDPAIAAPALAALLREEGREAAEAQIRHLHGRRQRRALVPLLPSAATTVSQPWRKGGVYLITGGAGGLGWLFAQEIAARCEGAGIILLGRSALAPGQATRIERLGREGTRVVYRQADVTDAAAVEAAVQAARTLGPLRGVLHAAGVLRDGLLLNKSEEQARAVLAPKLTGSLVLDRATRELDLDFFVLFSSISALGSVGQSDYAVANAFLDEFAHWRATAVARGDRRGRSLAIAWPLWAEGGMRPAVAALAGMPEAWRSVAMPSAQGLAYFYRALASVYAQVVVPSLAQAAPAEATSLVAGAVVASTVILPTAPAAPSSTPSTHKASMASAAVSMIDVEAFAARLRGVVLGLVTQLLGVPTEQVDLDEEFNAFGFDSISLTNLANRLNQELRLNLTPAIFFEHASVNRFVEHLLAEHADRLDFLVQTAAPARAEEAKSMPTVTSAVAATASDDIAIIGMSARLPMAVDLDEFWDNLLAGRDCIGEIPADRWDWREIYGDPLREPNKTDAKWGGFIDGVADFDPLFFGISPREAEAMDPHQRLLMSYVWQALEDAGYTRQALSGSDTGVFVGLGGSDYGQRVVAADGGVESHTLMGLLPSMAPARMSHFMNWHGPSEFIDTACSSSLVAVHRAAQAIAAGDCSLAVVGGVMAILSPTSHIGFSKAGMLSKDGRCKTFSSEADGYVRGEGVGMLVLKRLSAARADGDCIHAVIRASVVNHGGRANTLTSPNPSAQAELIESACRKAGVAVEQIGYIEAHGTGTRLGDPIEINGLKSAFKTLGHPMRGAWCGLGSVKTNIGHLELAAGIAGLLKLVLQMRHRTLVASLHCEQVNPYIELEGSPFYLVRENRPWPTGRDAFGRPAPRLAGVSSFGFGGVNAHVLVQEHLADAPPEEDEASVLIVLSARSENALRGRARDLLVYLERRGLDRSTASVYDNRRTELEQHIRLALAGLLDVDVHEISIDESFEAYGLDALARNRLAAALDESSPALLAATLRAGSVAHLADEWLASQGYAPVSVASVESSIRLGDLAYTLQVGREPMEHRLGFIAASFAQLNERLRAFLDGQDGMHELHRGQAKFGKNGLSLLGADDDMEGTIASWVVKKKYALMLKLWVNGMELDWRLLHVGRQLQRISLPAYPFDYRRCWVGETPAAARETDKAAGASEPPALDVALDETSGEVNQAESAQFLLPVWDAVHPELLADTPPQGRILLYGGDQGMVDAWRGGAVGLATLVVTPGASIEALSAQLAAAGSICEVLFVAPGEADVDPGNAQLIDAQQQGVLGLFRLVKALLATGHGLRDLAWTVVTTQAVDLGDGLPVRPAHAAIHGLAGALAREYPHWRLHLVDVAVDEADDWRNWRRLPPEENGATWCRRGGEWFRQSLLPYQPDTTNLPLPYRQGGIYLVIGGAGGLGEAWSRHMIERYQAQLIWLGRSAPNARIGERIDALAALGPVPRYIQADAGDEAALRRAAAEIAREFPRLDGIVHSALVLADQGLASMSEAEFSGALAAKVDTAVLSTRVFTNQPLDFVLFFSSMMSFGRAAGQSNYAAGCMFADAYARRLAQLARYPVKIVNWGYWGSVGIVSRADYQARMRRAGIGSIEPAEGMRALEVLLRAPVPQLGFVKTRIGVGGDSLRQYADTIPAVTARVMLDTPVPALRREDLLFADQRLESLVRGLLARQLETLGEHPVLPLYQRWLAESRRVLAEVPPPATTVEALWRAWDKHALTRRHLDGSTPELDLLRNCLQALPAILLGQRTATEVLFPQSSLSLVEAIYRDNPVADLFNGRLQQWVLAYLRERLARDPRAQLRILEIGAGTGGTSAGLLAALRPVHEHIAEYSFTDVSRAFLIQAQERFSVGFPSLVTRLFDVEQPLAAQDLPADYFDLVVAANVLHATSEIRTALRNAKAALRKGGVLLLNEIAGQSLFTHLSFGLLEGWWRYRDPAVRLPASPALSPKNWARVLREEGFGEIDFLLKERLDLGQQLVVAQSNGLVRQAEARGGGAEKGRHGSVVVSNTAFVATQARPAAPAVPPQAVAHATVEEPVAVDAVVHRVLLEELARTLKLDSGEIDTQMPFSDYGIDSILGVGLVKQFNDRLGLNLNTTILFEHSSLERLVVHILQRHRREAITALGLVVGPVLGAAEPAPPPAPAVAASCESGRASGAVDELIPAPEVELCRAETGPLEIAVIGMSAQVPGADDAETFWHNLMHGVDGVTQLPAAYLDRARQGDDKQAYYLWGGTLAARAAFDPQFFSISPREAKSMNPHQRLVLQESWKALEDAGYDPRSLAQRRVGSFVGAEPSGYFHETFTGSSDAIIASRLAYFLDLKGPALVINTGCSSSATAIYLACESLRHGESELALAGGVYAVLNETGLVSLAQLDMLSNSGRCHSFDAAADGTVFSEAVGMVVLKRMDDALRDGDPIYASIVASGINQDGASNGITAPSGEAQEQLLLETYRRFSIDPRHIGYVEAHGTGTRLGDPVEANALLRAFRAHTTDRDYCAVGSAKAHIGHTAAASGVIGLVKLLLSLKHGRTPGLLHFRELNPLIEWQDSPFFIPTANRDWPVESGRPRMAALNSFGHSGTNVHLVLKEYVDPVEAADDAQQPGWHLLPLSATNETRLRAYAGKLARYLREAGEGVNLGDVAHTLQHGRIALKRRWCLLVANAAGAIAQLQAYADGADVAELAGVHGLADGTLVTDRIEAAATPELLAALARRWASGAMVDWPLAGPTARRIHLPSYPFARDHYWMDGGALAATTVRTVTAAPRVPIMHPLLHARVDTATGPRFSAQLSGEEAFFRDHRILGRPTFPGAAYLEMACAAAARVLGTEALRLSAVVWSRPLQASAGGVALGIELQARDGDGCWRFEIAAAGQVHCQGLAAKAEARENAFERQDLGALLARMDRQGVGHDDCYAAFDAVGIHYGPSHRAIAHLHLGHDEVLARLVLDPAAALDSEATAAAYRLHPSLLDAALQATLGLSVGSARGGAYVPFALERLVLCASTGFQLWAWVRYTAGSRADEAVRKYDIDLLDDQGGLIARLLGFAFRPVEAVAPASADPVLLWCRDWKPDPVPLVGDSAPVRLSVLVSDRAAWRFDEPTDGSLLAHLQIDSIQEAPDQWYGRVLRGLFELIKAQAQKTVQACLLQVLVPAWGRAALLSGLSGLLRTATLEFPRLSTQLLAFDAPPVGLSALLEANARQPWHAQLRYLDGTRLAPAWAEVAGSTGAQMSALPWREHGVYLLTGGLGGLGRLLAEEILGQTPDATVVLCGRAAPDAEAAQWLRTTGGNGRLRYRRADVADAAQVREMVADIVVCCGTLHGVFHSAGVLRDSYLPGKTSAQIDAVLAPKVAGTVNLDEATQSLHLDLFVLFAAGAGALGNPGQADYAAANAFLDAYAGWRNEQAAQGRRRGQAVAFDWPLWAEGGMRISDSQREAIEAANGMRALERDEGIRALYRGIASGEAQLLVGVGDHARLRAWIQHLHAEPTAPAARSAQTGTVVPAGTVRALRERTRQFLIDSAASILELQPQDLDPRDELSDYGFDSITLTELTHRLNRQFELELVPTVLFEHPTIAKLTAHLLESFPEALARAFPEDVAGPPPRATQAPLATPERDMGEMRAADALRNTNTNTHAGAVAIIGMSGRFPGAPDLASFWQVLAEGRDCIGEIPADRWDWRQYYGDATQDGRHTRVKEAGFIEGVAEFDPLFFGITPREAGLMDPQQRLLLTHAWAAIEDAGYAPSSLAGGRTAIFVGTAPSGYSSLIEQAAMGLDGHSSTGSVASLGPNRLSYLLDLHGPSEPVETACSSALVAIHRAVRAIRHGDCETALAGGVNTIVLPEVHISFDKAGMLSPDGRCKTFSRHADGYGRGEGVGILLLKDLAAAERDGDAIHAVIRSTVENHGGRASSLTAPNSKAQTELIKAAVRDAGIEADSIGYIETHGTGTALGDPIEINGLKAAFEELQAEAGGAAPIPDSCALGSVKTNIGHLELAAGVAGVIKVVLQLRHRTLARSLHAEEPNPYIRLAGTPFYLLGETRPWPMPRDAGGRPLPRRAGVSSFGFGGVNAHVLLEEYVGSRAELPASEPRGPYLFVLSARQPARLVEQAAQLVAFLRDETRPALAELIYTLQVGRDAMEERAAFLVADYAELTQRLAEFVAGRAGAWLHRGRARHDSVAATLFDDQNEQREAVRVWLAQGRHDKLLGLWTQGLEVDWRVLHAGPPPRRAHLPTYAFARERYWPPRPDRPAARPETRVNPLMMQEEAANEDELDRHSALLDQLIHRQISADEAVRLARRQEQL",
        ],
        "Ia6RrYNflQpEjxBCKTb5azk9_FTDvB-5": [
            "AXA20092.1",
            "AXA20092.1",
            set(),
            "MTDHVAEIYHQVASGQLSKDEALARLCQVKDEQAAQGVAKQVHKPAQPVAVSREAVTQALAALYAGQSKIAPEAIDPLEPLESYGIDSIVIAGMNRELGERFASLSKTLFFEHRTLDSLAGFLWREREPACRAWLLASGTTVRTVAASPASRIEVAASAAAASAISDEPIPGDRAGTGAPETFGPANGMAVDASAWPPIAIVGLAGRYPQAADLDAFWRNLREGRDCIVEVPVERWPLDGFYEPDPDRAIASGRSYGKWGGFLEGFADFDALFFSISPLEAMGMDPQERLFLQSAWHALEDAGYTRESLAERCGGRVGVFAGVTRNGFGLLGLEAWHAGGMVFSQPAFSSIPNRVSYALDLRGPSLPVDTMCSSSLTAIHEACAALHRGDCVMALAGGVNLCVHPASYVGLATGRMLSRDGRCRSFGAGGDGYVPGEGVGVVMLKPLAEAQASGDRIHGVIRATSVNHGGRTNGFTVPNPTAQAELIAESLRKSGIHPRAIGYVEAHGTGTALGDPIEVSGLVQAFAPYTRERGFCALGSAKSNLGHLEAAAGIAGLTKVLLQMRHGELAPSLHAAQLNPNIDFEGTPFVVQRELAPWAEPELDLDGQLRRYPRIASVSSFGAGGANAHLLVEQYLAPAVPSVSSAGPQAIVLSARTPERLCAAVQALLDHVESGGGTAAGLAANLSAWLVEELAAIVGVEATAIDPHETLANLGVEILHRTRWYERVQERLNLPWSLKNFLDQDSVQQLGNTLLREQGATVVAQFHAPVAAPVLADLAYTLQIGREAMPERLAFVAADLAELATGLRGFLDGASRLPLWHGKAARGRALPARVDQAQWQAWIDARDWSQLLPAWVAGHELPWRDMPGAPGARRIGLPLYPFAAERYWVDPQSLRPRPTASLESRLHPLVRKRVDSVEGPAFLSRFSGAESFFSDHRVGGRSILPGVAYLEMARAAASLAADGATIRSLRNVVWARPIEAGADGVAVTLRLQPHQQGSWRYEVLGADDGVHGQGLAELALPEAPPVDLDLAALRARMQGGALANEALYQAYAAMGIAYGAAHRGFVQALVGESELLAELCLPSAVQADAQAYVLHPSLMDAAFQATLGLYLLSKRADAAKAMLPFALETLELHWAPPARVWAWIRSRGERSGIEKFDIELCDEQGRICVRMLGFSSRVLEAPAVSPEVPPAVLEAPSLLLSRYAWQDAPARRAAPDPALTRRLLLVSIAPQPVVWRQLGQGELLQAAAVQPEQACASLYVQIFERVQAWLEEKGRDTCLWQLAISGQGAELLLAGLSGLLRSASQESRRLLGQLMILEGDEDLASLRARLDENAASPFDSLVRYRAGRRETWHLLELPSNDGEAEAAPLPWRQNGVYLLSGGAGELGLLFVEEIARRATGATLVLTGRSALPDARRARLDALCEQGAQYRYEPVDVTDRQAVTQLVEYVVAEYGRLDGVLHIAGVLRDSYILKKDRAAFEQVLAPKLLGTANLDHATRTLDLDFFLMFSSSAAIFGNLGQTDYAAANGFMDAYAAYRQARGGRGRSLSVNWPLWRDGGMGMEAATEEMMLANTGMVAMRTPSGFSALARALHSDLPQVAVMEGLVERMRQKLLVPSAPSMQAVPVASASAAIAPSAQTDHHAEIVARVARGLRQMVAELLKLELDQIDIEDDLSDYGFDSITFTSFSNRINKQFGLELIPTIFFEYPDIAGLAGHLAEAHGAALGASLGLLAASAQGDRAQTRSAMSAEAVATANVSAEMPVPLAPQLSDQSAAASNTPVARRGVAVIGISGSFPGADGVDALWQLLERGGDAICEVPASRWDWRRCLPPGESEAVQARVRWGGFIDGVDRFDPLFFGISPREAELMDPQQRLLLSYAWLAVEDAGYAPQSLGGTDTGLFVGTAVGSYGSLVVQAGRSRDAYSSTSSVASIGPSRVSYFLDWHGPSEPIETACSSSLVAVHRAVQAIESGRCEAVLAGGVNTISTPEAHIAFSKAGMLSVDGRCQTFSAKANGYVRGEGVGMLFLKDLVAAERDGDTIHAVIVGSAENHGGRATSLTAPNPKAQAALLKAAYAKAGFDPRLLGYIEVHGTGTELGDPIEVNALKSAFKDLYQRAGVEPPGQPHCGLGSIKTNIGHLELAAGVAGIIKVLLQLRHRTLVRSLHGEQVNPYVQLEGSPFYLVQENLPWDAPRDAQGREQPRRAGVSSFGFGGVNAHVVLEEYVAPPARATAPAACPVLVPLSARNETRLREAAARLADFAAAHADDAALDLHDLAYTLQVGRDAMEARLGLMVSDKAELARCLRAWLDDAGAGEVFQAAPGKAQKEALALFAGDEELAGVVEGWWRNGKQAKLLDLWVKGLDLDWARLRAGAGRRRISLPGYPFANERYWLKPVSAETAALGLVASTPIETDEAAVLFFEENWHPYPIREALIASSTRTLLCCLSDATHRKALREAVFRYDPKWHLVFLDRQAPEPFDRQGWTDALCLLEQGGPVIDAVLYLRPLEEAGLRLAEAAPLGLVQALGSMKTRPARLVLGGEYADESERSQLEAWIGLERSIGLALPGCRAVTVLREAGDTIDWVAWTQLLRTALGEAAPRNLLADRDSLRHLQVQPLEPRAAAVADLGTTVLITGGTGGLGLILARHLAVGRRCNLVLVGRSPFDVVRQAAVQALQAAGSEVLYLSADVADAVAMREVVAQARARFGSIDSVIHAAGIQHAVPLADKQSEDMRRVLDPKVRGALVLDQVLAGEPLRLVCYFSSSSSVLGDFGSADYALANRFLSAHALARERRRARGERAGRSLSIEWPLWREGGMGVGDDAGTALYLKSSGQRLLEQTEGLAAFERLLASGATRALVLVGERERLHRMLGLAQVPSAPATQAMAVLMPSQTQTFSTASLEEQVSAELSVLIGDQVKLAPELLDAESNLADFGLDSFGLAELARALSARYDIEVAPSIFFAYSSIARLVGYLLDKHRAEVQAHRHRTATRVDVAAIAPQPASLAPPAAISMSTPAPPQLPVANAIEPAPTVPAFDGEREPIAIIGISGRFPKARDVDQMWRILAEGIDAVDEIPVERFDWRDYYCGLEAQPGRTNSKWAGCLDGVDEFDPLFFEISPREALAMDPRQRLLLQESWNALEDAGYGPHQLRAGPVGIFVGVEEGDYQRVVPDPGVTSNHNGILASRLAYFLDLNGPVLAINTACSSGLVALHQACASLLSGESDTAIAAGANLILTPEPYIGMSQAGMLSPDGRCRAFDRSANGMVPGEAVAVVVLKRWSRAVADGDPIRGLIRASGINHDGRSNGITAPNALAQASLVRQVQRAAGVLPEQIDYVVTHGTGTRLGDPVEIQALVEAFGQPVDGRAYCALTSSKGNFGHTFAASGLLSLIGLVKALEYDTIPPSLYCDQDSDYIAWRDSAFRVNKQARSWPRASGRARLGAVSAFGMSGTNAHVLVQEAPVAVARVAVAQTDVVLALSGKTEAALRERLFGLRDWLASAAAERCELASVSRTLLDGRHHFAHRAAVVVANRASAIAALEHLATLDTADGPDYYLGKAPRGFKGEVALRERAANWPALPSVDAAAYRERLGELARLHCQGYTVAWSALDGLQPAARVHLPGYPFARESYWPKTRISPPVATSVPASPAFPSSHSTPVVISLRPMLRTELSDQARGHARACYVARFSGEEFFLRDHRVRGQAVLPGVAYLELARAALEASSGRPVPSGLQLRHVTWVQPLMVDEPGVEAYIDLHRQDNGEWRYELASGAHDESERYLHGQGFLALVAQAEPTALDLSVLHGNCRVTEFDSAACYAAYQAVGIEYGPSFRAVQRIWVGEGQALVQLRLPVEALSDSGAYTLHPSLLDGALQASIGLAMAQATQGGEPMLSLPFALDSLVLHWPCPNETWVWIRPTPGAQSSRVRKLDLDLCDAQGRVCVALRGFSARLVAQGGTAQPALIPARVEAERVSMAAPINGLASASLPAKAKGVVPILAPAAKATGASTALPATTEGNASSAIPVVKTGAAPGAWLPVGLTMLAPLWTVRRVDDGSEPPAPVHMLLIGGDATQRAIWRQAYPQLRAIDVAPSTTIDELRGQLAATGVIDELVWIAPAQHSQDPTDEALLTAQASGTLALFRLVKALLAEDYASRRLAVTVLTRATQQVHPQDLVAPAHATVHGLAGSLAKEYPHWSVRLIDLDGQNADPPPERCRALPVGSESWAWRRGEWHKPELIALDEPTRGKVPAPYREGGIYVVIGGAGGLGEVWTRHLIEHYRAQVVWIGRRAEDAGLRARLAALAAHGSAPVYLQADAGNRLALSRARETILQRHGRIDGVVHSALVLQDRSLARMDETTLQSALRPKLDVSLRIAQVFADAALDFVLFFSSMMSFSRAAGQGNYAAGCTFKDALALALARRWPGAVKVMNWGYWGSVGVVADARYQERMSRAGIGSIEAEEAMVVLERLLGGPDAQLGLIKLSRAQAVEGVRDDLRGARYGAALPALLPQLAARPLPPEHASRLAAAQAALPPQAMQALSLRLLGQALLGFSLDGRHLLPGGAAGLALAGHYRRWYDTSLRLLDAGGWLQSLPNGDYQILATAGQQDAWPEWEAARAAWLANPQQQAWAQLLEVCLRALPELLTGQRKATDVMFPNSSLRLVEGIYRGNPIADLHNHILFDALEAYVLERLAREPGTRLRLLEIGAGTGGTSAGLLQRLDRYAANIDEYCYTDLSKAFLLHAEQHYAPGRPFLRTKRFNVEEPPQAQGIAADSYDIVVAANVLHATVNIRRTLRHAKVPLRAGGLLALNELGELSLLTHLSFGLLDGWWLYEDPALRLEGSPGLSSEGWERVLAEEGYAPLWRPAEECNRYGQQVLLAQSDGRVKRDVAVPPEMAAAPVEEVSAPASAFTSAQVADSTPAATFAPVTAPIPVPDSRDLEQAVADHVRTLLRECIGKGLDLDPRRIEADRSFSEYGVDSILAVQLVNEINQRLGIVLQTTVLFDYSHLDVLAEYLEQTHQAALRASLPEVSEVPALQAASTLAPKIQAQPSGVAFPLISPTQPIGATLPFVPPSVTGSHRRALISGPGQIQDLRLVAMEVPASLQPRQVRVAVFASSLNFSDLLCVMGLYPNMPAYPFTPGIEASGLVLEVGSAVSTLCPGDEVVCLAQGCHATEIVCHETQAWAKSPQLSFEQACALPVVALTMIDAFHKADLQPGECILIQTAAGGTGLIAMQLARHYGATILATAGSQEKLDYLRDQGAQHLINYREQDFEAEVARITGGRGVDVVINTLSGEAIDKGLRSLSPGGRYIEIAMMALKSAQAVDLSVLDSNQSFFSIDLARLIAERPEKLEQYRRELASLVEQGVLLPTMSRVFALDQLHDAYRYLQDRRNIGKVVLQVPQAVPLADQASAARAVDAVAGVHKAQPVPVNYADEPIAVIGMSGRFAHSPDLDSFWSHLAKGHDLVDPVLRWDLSPSGGRCRDGSFLDEIDRFDPLFFRMSGLEATYMDPQQRLFLEEAWHTLEDAGYAGEAVKGKLCGVYVGCTRGDYAQLCKSAPPQAFWGNSGALIPARIAYYLDLQGPAVAVDTACSSSLVAVHLACQGLRSGDTELALAGGVFVQSTPGFYLAANPAGMLSATGRCHAFDESADGFVPGEGVGAILLKRLSDAIADGDHIHGVIRGGAINQDGRSNGITAPSARSQERLERQVYDRYAVHPETLQMIEAHGTGTQLGDPIEYRALRQAFGHYTQRVGFCALGSVKTNIGHLANAAGIAGILKILLALRHRQLPPSLHFRKGNPAIDFEGSPFYVNTELRLWPAGERAPRRGAVSSFGFSGTNAHLVIEEAPAVPLVARATRRELELVVLSARTAGQLREQAARLLAHCQAQPQTSLGDLAYTLLCGREHRGYRLAAVVRDLAELCEVLSAWLEQGDDSRLQLGALDESGVREQLQQRRLGQVAIETVRAGQLEKLSSVAELFAQGYKLDYAGLFGSGYRRLALPTYPFAQGRYWVDDSLQHAVVPSTPATAPVVPTPVVPAPAVTQAEARQAPSRISTVPDQVMREASVAYLKQLVATTLRVSPTEISAHEPLERYGIDSILVVQLTDSLRQHFDSVGSTLLFEVQTIDALAERLLATEAPALARQLGMDAVAPLEATGSAIEAELPESPPPQPETNSQVQPAPAVVTVAETVGASAAAESSQPKAHGDVAVIGMSGRYPKAIDLNEFWWNLRAGRDCIDEVPAQRWDWRKHFDAQRGLHGRSYSRWGGFIDGVDQFDPRFFRIPPSEAEHIDPQERLFLQTAWLAIEDAGYTPSTLSAKRRVGVFVGAANSTYTLLPSHWSIANRVSFALDFHGPSLAVNSACSSSLTALHLALDSLAHGSSEVAVVGGVSLVLHPMHFNRLSSLGMLSSDAHSRPLGEHADGFVDGEGVGALLLKPLQRAIEDGDSIHGVIKGSMVNAAGKTRSFAVPDAAAQARLVREAQARAGVEADTIGYLEAHSNGGELGDITEMQGLAEAFAGTAERGHRCAIGSVKSNIGYCESAAGIAGLTKVLLQLRYGELVPTLHARCANPRIDFAGTPFALQQELSAWPRPANHPRRAGVSAFGVGGAYAHVIVEEYVAPVETQPEATGRALPIVLSAANAERLRVLAKRLAGFLGSEAGRRTALTDLAYTLQVGREPLAERLGFIAESVEQVREVLLAVAEDREVPLPLVRASLDRGRAGWAMFAEDEDFKRTVEQWIAREKHASLLDLWCRGYPLDWRHLYAAHRPRRIGLLPGYPFAEESYWAPESLRYAGVLEDADAFDATPFEPDQPAGEQS",
        ],
        "RyDIaUZc_b21_kQalx7J3yNO4l5f-439": [
            "AXA20093.1",
            "AXA20093.1",
            set(),
            "MKPNLNQDFDATPSSHAADRRPQAGAGGCVRAPGHVDTAIIGISARYPKAVDWRQFWENLRAGRDCIVEIPPERWDWRAYHDSARGTPGRSYTRWGGFLDGIDRFDPRFFRIAPSEAEHLDPQERLFLESAYLLIEDAGYTPASLSASRRVAVFVGAMNSSYSLLASQWTLANRVSHVFDFHGPSLVVNSACSSSLSAIHLAIESLATGTSEVAIAGGVNLIMHPAHYARLASVGMLSAGSHCRAFGAGADGFVDGEGVGSVLLKPLQRAIEDGDLIYGVIKGSALNAGGRTHGYTVPSPVAQGRLVAEAIERAGFAPHSIGCVEAHGTGTELGDPIEVRGLAEAFGAPVGAAPWCALGSVKSNIGHGEGVAGIAGLTKLLLQMRHGQLAPSLHADTLNPRIDFNGTPFSVQRKLAPWPRPAGHPRRAGVSSFGAGGANAHLLVEEYVAPIVPAPTDADSPALIVLSAANLDRLRAVAQRLLDFLNGEFSSGITLAELAYTLQVGREALAERLGFVADSLGQVRACLAAFLEDREAGRPLLRSSVGQGRAVGSGMLDDESFAQTLRGWIKRGKHELLLKLWGQGEPLDWSLLYRGARPRRVSLPGYPFAGERYWAPAAVRYAGVVCSRRRPAIDPDLLLCQPTWRAASLPAATGRALPHRELWLLGSQARLDDAALPALPIERFRSEQAEPVARCVDLYGQFHARLRARLRDKLSEPLLLQVAILGRDDELLLSGLSSMLRSLGQESRKLSGQLLVMEGGEDRATWQARLDENAVRAHEDWVRYRQGRRETWALQELPPATVEPALPWRARGVYLLSGGAGGLGLLFAEEITRRAEGATVILASRSAPVETRRARLAALAEQGLAIRHAVLDITDAAAVQALVDEIVASYGRLDGVLHLAGVLRDAYLVDQERDRVDQVLAPKLLGALHLDLATRMLPLDCFVLFSSAAALFGNAGQADYASANGFLDAFAVYRQTQGRSGRSLSVNWPLWAEGGMSMDKATEQLLTAGTGIRPMRAATGCRALAHALAGAFPQVLALEGEPVHMRAALLGQSSPVAVATAGDATSPKTSLSCQVRELVAALLKVEPEQIEAGQDIGDYGLDSIGFTHLANQLNLRFGSALRSTDFMELEMASVERIARLLEQKLPALSTGTTPVTRRDRGPSAVPVTPTPRDAGPAAHSDEPLRAQVREAVAALLKVEIEEVDLDLDISDYGVDSIGFTHLANRLNEQQGTRLRATDLLELEQVSVIRIARLLREDPGSRTLLDAVGADASLAVVEGR",
        ],
        "DTee9G4M8sEfnM4HaPfI37rT74pq7M_G": [
            "AXA20094.1",
            "AXA20094.1",
            set(),
            "MSVYMFPGQGSQAIGMGTDLFVSFPELTEAADRILGYSIRELCLEDPKHQLGQTRYTQPALYVVNALSYRRRLHEYGAPAYVLGHSLGEYNALEAAGVIGFEDGLRLVRKRGELMSEAPPGAMAAVIGPDEVAISALLARHGLDAIDIANLNSPSQTIISGLKEDIARAAPLFDAEQAHFVPLNTSGAFHSRYMTVARQAFVAYLGEFHFNRPRIPVISNVEAQPYVLERTAELLAAQITQPVQWTRSVHYLLALGQSEFLELGPGQVLTRLLVEIRKHTPTVAPATAGSLSSTPAYDRERELSELQQRIDDWNGRYPVGTRVQVERYPQQLVTRTPAMSLFGHRAAIYLEGYNGYFDLADVHPLHGASA",
        ],
        "mB22-i4RqtslyO7_HappM4rJ4Z2Qbkfn": [
            "AXA20095.1",
            "AXA20095.1",
            set(),
            "MNAIPKDYAVPSISIERLGSVAFKQDHGLRYPYVAGSMVKGIASTAMVINMGRAGFLGYFGTGALDAVSIERAILEIQAALGDRQPYGMNLLSNASTPQAEMDTVDLFLKHGVRRVEASAYMQITVPLVKYRASGLRRDAQGAVVARNMILAKLSRPEVAALFLSPPPDKLLAELVAGGAISTAEADLARLLPMADDICVEADSGGHTDMGVLSALLPSIVRLRDELVAHHGYARTVRVGGAGGIGTPEAAATAFILGADFILTGSINQCTVEAGTSEAAKDLLQQVNVQDMDYCPSGSLFELGAKTQVLKKGVFFPARANKLYELWKNHSSWEEIDAKTREQIQNKYFMRSFESVYEETRAYFLRAEPGEIEKAEKTPKHKLALVFRWYFVHTMRLAMSGSSQQKVDYQIHCGPAMGAFNQWVKGTELESWRNRRVAEIAHRLLEETVQLLNRSFLAMSS",
        ],
        "iI7aI2dI9vaha9f0rVTi_YFrfMXjY1eh": [
            "AXA20096.1",
            "AXA20096.1",
            set(),
            "MLELTKRLADALVSISLFAACRESGLGALLKDRAGLPTRAAQLTWLAPQCGIDEARLEAALQALRDAGWIEALEDGRLIPRATFERVEPWSEAVAVGLDRDWGALLREQDGRRLRHWLEQGAAARESLAGCQAEAEALDAAAMAPLLFELARLDDAAWLQGRDVTSLAPANAALLRADFLRRGWSLDEAAGLIPNVQGLAMLRDAAALGPLLFLARRPEAGARTLASTVALYRANLRWRNALDQAIAVADSSAHEAWPTGACVMAAAAIGCLPERPSPSRPAVKPGPARFALHQWTARPYRVRHPSLDDLAILRELDLASWPVGMAVPENELRRRIEQFPQGQLLIEQDGEVIASLYAQRIDTLDQLRHTPYARFAWIHRPRGALAHLMGICVAPDWQGHGLADQLIDFCLVYLASFEGIDSVAAVTRCHEYGRFGDKVTLDDYIRQRDEEGRYREPMLQFHASHGAIIHEVVPGFRPEDQANHGTGVLVEYAHYRQAPELAGPVVAVDSVGPGSAALDVAEAVRASILDVLGELHAAAYGPQVPLMEMGFSSFHLQELQRSLGERVGLKLDATFFFQHGTPAGIVEHLRERLAPIGTEQADIRSTVDTETTPSSATTDGIDVPERIAVIGVACRFPGGVGNPEQFWTLLENGVDAIGEREPGVSPTAASTRRGGFISAVDRFDAGFFRISPREAELVDPQHRLLLEVVWEALEQASIAPGRLAGSDTGVFVGVMGHDYERLLRQQGGAPPIDPYFATGNANSIAAGRIAYYYDWHGPTLAVDTACSSSLVATHLACESLLAGECSLALAGGVNLLLHEDMFAAFEQAGMLSPEARCKTFDASADGYVRGEGCAMLVLKRFSEAQRDGDPVWGVIRGSAINQDGASAGLTAPNQGAQQAVIEAALRKGGVVPHALRYLEAHGTGTRLGDPIEVLAACAALSAGRPIGQPLLLGSVKTNIGHLEAAAGMAGLIKVMLSMRHGLIPRHLHLQQPNPHLDWAALPVEVVSEARPWPVGPKLAGVSSFGFSGTNAHVVLEEYPANPANTVPMAARSSALLLLSAKREEVLQTQVRQLHEAIGALDEADLPDVAYTLQVGRDAMEYRLALAVGSLAELRQALARFLAGEAGIRQLWQGRAGQQGYLLGSFVLDEAFTASIATSLAAWFARGELGKLAELWVQGLDVDWRRLYGANPPRRISLPTYPFMKERHWLPQAVAQATAEASGAPLLHPLVHRNTSNLAEQRYSSRLDQQAFYLRDHVVQGRHVLPGVAQLEWARAALALALGDTSASLRLEQVSWVQALTVEQALEVHIGIEADEGGWLTYEIYRGSDDEVELYSLGRARLDAERKVPNLDLATLQARCTRRIDGPACYARFTRMGLGYGPAFQVLTELHVGADLAIGRLQVPVGIELGDYRWSPSLLDGALQASYGLVDETAGLQLPFAVESVEQSHALPESALVVVQRAADDSGVLRKLDISIVEESGRVALRLTGFSTRAVQAAAPADSLLMVPRWQARSAVEAPPEPGYRTHRVVLCEFEALRGGFDAALPAASVVHWQAPGSLAERYARYAGQLLCELQALAADHPADPVLLQLVVPAQGEAAVLQGLVGLLCTAQQEYPWLHSQVIALPADAPVADCLAREAAAPVPRVRYQGMQASTREVMDFVEVPPTETAMPWREDGVYWITGGLGGLGLLFAAHIARQVQTPVLVLSGRREPDTAGQAQLDALRALGANVEYHALDIADAAAVAALARNIIARHGCLNGVIHGAGVLRDGLLHSKTVDELQEVLAPKVTGMMALDHATADLELDWLLLCSSMTTVIGNTGQGDYGAANAYLDAYAVHHEQLVAQGLRRGRAISVSWPLWAEGGMRIDAEGQAYLRRSTGMQSLPSEIGLAALEQLLATPRAHSLLLYGDRTRLLARVQALYRAPEPVVVRSVLALAPAAGPVDSQEALRKTARRYLTRLLSRSLKLPPQRIDVQTPLEQYGINSILVVSLTRDLEASFGRLPATLFFEYQSIAALTEYFLAHHADSLTMLGAAPAATQLMAVTSSAPAREASIDAPALRRRRHRRSLPGSMVAGPPVSTVGAPLDIAIVGLSGRYPQARSVADYWANLLKGIDCVTEIPAERWDWRQHFDARKGQDGKSYSKWGGFIDGMDAFDPMFFGIAPREAQLMDPQERLFLQCAYHAIEDAGYTRAGLAASATEGERRGQVGVFVGVMYEEYQLYASQAQARGQGLSLFGSASSIANRVSYHCNFHGPSLAVDTMCSSSLTAIHLACQSLRQGGCSVAIAGGVNVSVHPNKYLMLSDRQFMASNGRCTSFGEGGDGYVPAEGVGAVVLKRLEKAIADGDHIHGVIKGSALNHGGKTNGYTVPNPVAQGQVISQALAESGVPARAISYVEAHGTGTLLGDPIEIAGLSQAYGASTQDKAYCAIGSAKSNIGHAESAAGMAGLTKVLLQMQHGQLVKSLHSDTLNPHIDFSQTPFVVQRELGPWTRPVLEVDGGGEHEYPRIAGLSSFGAGGANAHLIVSEAPASSRHEAVPRQGPVLVVLSARNENILRCQAEQLLTHVQSHASDLANLAYTLQVGREAMEHRLAIVAASTEQLSARLNAYLQDDTLDEAVYRGEPRRSQEAMAVFGGDEELQEVVAKWIARGKLEKLAELWVQGLLIHWEQLYGQAMPRRISLPTYPFAQERYWIDAGAATMRVAGAQILHPLVHVNSSDLQGQRYTTVLDAGTGLLRGHRLHSRPTLPALAQLEWARAALAHALGGTAGLCLEEMRWLVPLHVDAPTTLHIALDWEDETHVGYEIYREDDEGREVYAEGRAELVDALPAPRLDLPALQAQCTQHLDGDEAYSRLATAGWSCADSFQALSSLQSGEGLAIAHWRQKVDASWQDYALAPNLLDVALQACRLAWPQQDWSWPTAARQLRMVGSLPVQGMVVVRQHPGQLDVDIDFADGEGYVLASLQGLAPQQPSTQASPVAQTLLLAPCWTPQAGPPAACERPSYAAHWVVLCELDAPASLEAELVPAHCLRWQAEGNPAERYGVYAGQALAWLQKIVAGGPSGQVLLQLVLPARGEAALMQGLGGMLRSARLEYPWLLIQVIAVDSAQELAVRLNTEAIAPVPALRDGAAGREVLDFIPLAPPEGVRAMPLAWRDEGVYWITGGLGGLGRLFATAIAAEVRCPVLVLSARRPPDTAQAAFIERLREQGARVEFRAMDTGDAAAVEAVARAIVAEHDGLNGIIHSAGVLRDGLLANKREADLRQVLNAKVGGLFALDMATRDIGLDWFLLCSSVSSVLGNAGQTDYAAANGFMDAYSAYRQELVEQGRRRGRCVSLSWPLWAEGGMHIDAVAQEQMRRATGMQALPRAAGLAALHQALAATVSHVLVLHGEPQRLRDYVSTAYRVPALAEPVAKMSPGRGYRRELKGLDLADCVNWDLVEHTSALLQMPRDAVDTQANLIDYGYDSVSLTAFAARLGEHYGIALTPSLFFSHPTLEQFSVYLLDSHGDALAAFYRAASEVERAEAGPAQLGAPTAATSASIRRRRHAQLAGAVANIVEPIAIIGISGRFPGARSVDELWTILRDGREVLQSAGTERFAAWPPPQRPACDRIGLLPGVAEFDPLFFEISPREAEAMDPRQRLLLQEAWRALEDAGYGVTQLQLHTVGMFVGVEQGDYQLIGKTEADVTSNHDGVLATRLAYALNLHGPAMAINTACSSGLVAAHQACLSLRVGECDTAIAAGVNLLLTPAIVRSMEQAGMLSPEGRCHAFDRRANGMVPGEAVVALVLKRLSRAEADGDPIHAVIVGSGINYDGKTNGITAPSGAAQTRLLQSVYARHHIDPADIDYIVTHGTGTPLGDPVEINALADAFTPHERAPQSCALTSSKTNLGHTQAASGLVSLVGLVQAIRHETIPASLHCEQLSDHIAWQKSPFYVNTAARPWPAPTTRARLGAVSAFGISGTNAHMVLRGHAAPADAGRHAAVRPLLLAVSAKTAEALRQRVQDLIERLQAREHDAAELASISHTLLVGRQHFAHRCAVVVQDREDAVYALQQALSRETRANLFRGVVSREFAAQKALLDYGQELIGRIAGVQQTTKQVQQQAGAQGEASREALSALADMYCQGYALAWHELFGQTPPNRVHLPTYPFARETYWVKPTKHAEAGEAVQLHPLVHRNTSDLDEQRYSSRLVVDAFFLRDHVVRGCSVLPGVAQLEWARAAVALALGGEPSIRLGQVDWLQPLVVEQAAECHIALAPLDDGRLAFEIYGDNGQVHSQGWAEAVSPGQVPRIDLAGLRARCTYRLTGEQCYARFVRMGLNYGPSFQSLAGLRRGEGIAIGELRWPADVDQEAAFVLPPSLLDGALQSCIGLYAESTGLILPFAVETVEQWGAVPATAYAVVQPGADDNEAVRKLDIRIVDEQGQVAVCLSGLSLRSVAPASTAVGTLMLAPRWRVQPALATNMVPANVAHCVIFCEVAPVDLRETLPTASSMHWTAEGSLDERYTRYAEKLCIELQTLEASRSDRLWLQLVVPAQGEHAVLQGLDGLLRTAGQEYPWLVAQTIAVKDTSNLAARLAAEASSPAPRLRYGEAGREILDYAEVFEPRQGVRPWRDRGVYWITGGLGGLGRLFAAHIARQAQVPVLVLSGRREPDAAGQAQLDALRALGAHVEYHALDITDVAVVAALAQNIVGRHGCLNGVIHGAGVLRDGLLRGKTVDALQQVLAPKVAGMMALDQATATIELDWLLLCASAAGVLGNVGQGDYAAANAYLDAYAVYRDELVAQGHRHGRAISVSWPLWAEGGMQVDAAMQAHLQRSTGMQALPSEAGLAALDQALSESGAGQVLVLHGERARLLGHVQAVHTALALEAQEETATLVAQEAGMDFREAAQRFLTRLFSQSLRLPPQRIDAKVPLEQYGINSILVISLTRDLEASFGRLPATLFFEYQSIAALTGYFLEHHGAALHALLGWSKAEAGAVSTLPTPSPSPAPVNAAPSLPARRFSGRLLDRFSRRNALTSGAPPDTPPVPLDIAIVGLSGRYPQARSVADYWANLLQGIDCVTEIPDERWDWRGQYDPAKGKLGKIYSKWGGFIDEVDAFDPLFFNISPREAELIDPQERLFLQCAYHAIEDAGYTRAGLAASATSGERRGQVGVFVGVMYEEYQLYGAQAQAQGQALSLFGSASSIANRVSYHCNFHGPSLAVDTMCSSSLTAIHLACQSLRQAGCAVAIAGGVNVSVHPNKYLLLSDRQFMASNGRCTSFGEGGDGYVPAEGVGAVVLKPLEKALADGDHIHGVIKGSALNHGGKTNGYTVPNPGAQGQVISQALTEAGVPARAISYVEAHGTGTSLGDPIEIAGLSQAYGASTQDKAYCAIGSAKSNIGHAESAAGMAGLTKVLLQMQHGQLVKSLHSDTLNPHIDFSQTPFVVQRELGPWTRPVLKFDGGGEREYPRIAGLSSFGAGGANAHLIVSEAPASSQHEAVPRKGPVLVVLSARNENILRHQAEQLLAHVQSYAPNLENLAYTLQVGREAMEHRLAIVAASIEQLSARLNACLQEDTLAEAVYRGEPRRSQEAMAVFGDDEELQEAVAKWITRGKLEKLAELWVQGLSIHWDRLYGQALPRRISLPTYPFARDRYWVPKNLPSIDAASTQAAVLHPLVHRNTSHLGGLRFSTRLDPQSWLLREHQVQDHGLLPGAAQLEWARAALSLALEGAKVRLRQVTWLRPLLAEGEAELHIVLRVEDDGRISFRIYREQDGDTLVYSHGWAEALGEQTPAPTLDLAGLLEGCTRHWSREEGYARLEAMGLHYGKNFQVLMSWQIGDAAVVAELRAPDVERLAGYGLPPNLLDGALQASLGLAGEQVGLSLPFAVELVEQWGPVPSPAYVVVHRAVGDSAVVPKLDIDIVDAQGQVAVRLGQFSRRSVEALADAESSALVAAADAVREWTLAPAWDIADLDAHQNANQDGFQKQGTLVLGEADWIGSSGLRNLDWEPEASRECLAERLGEQGELQQLIWAVPSAEPHAALMGLRLVQALLALDYGTRSLQLVVVTRQAQAVWPQETADPRQASVHGLVGSLAKEYPHWRISLIDLPAQITQDGHQWLAQAAQAADSRGDARAWRDCRWYRQQLVPCRLPAVQASAYRQGGVYVILGGAGGVGVAFSEYLVRHYQAQVVWLGRRAEDAVIAEQRARLGGFGPAPWYLRADATDRSALERAYARIRQRHGTIHGLVHAAIVLADRSLAGMPEAVFAAALDAKAATTENFDAVFGTEPLDFQLFFSSLQSYTKSAGQSNYAAGCCHADAYAHGLRQRRPYPVKLMNWGYWGSIGIVSAEGYRERMAQAGVASIEPPEAMIALERLVAAPLHQVAFVKTTTAQMPPLLAFDPQARIELAKASPALSLPAPVALPKPDAAYAEQAAMEARLARLLWAKLTAWGWHGAQAPGLVPAYALWHQASLRLLGEQAQPVVTAADEAILSAQWQAYVHALRDDKALGAHVRLADVALQALPAILRGERAATEVLFPQASLNLVEGVYRNNPVADYFNAVLGERLQAHVQARLAQDPRARLRILEIGAGTGGTSEGLFRCLAPHSERIAEYAYTDLSAAFLRHAEQQYASLAPYLHTQRLDIERSPLAQGFEAGSYDLVVATNVLHATRDIRRTLRHAKSLLRAEGRLLLNEIEGTSLFAHLTFGLLDGWWLARDPALRIPGTPALSWNSWREVLMGEGFRPVLAPAFEAHRFGQQIVEAGSDGVIRVQAEVANAPVVASEVAVAAPSASRPATMPTRRSHAAVNAAPVAAVATATGGARGRQARVRQAIRESVLEALKMNAAQLQDDQAFMTYGVDSITGVALVNTINTRLGLRLPTTTLFDYSTIEQLSTHISMQYAVQLSDAELEPVAAVVSAPVMTEMAASPLSESALDTVPMPSLAAEVFRAAAPVPVWRPSPAEPVQVQPRLLPIVAPGPSGSGPTYHRVWLDRPGSIDEVRIVADSLSPLQPHEVRIAVHAFSLNFGDLLCIKGMYPTQPAYPFTPGFEASGEVVAIGVGVSSVAVGDAVMAIAGAELGAHATVLTCMEQQVYAMPRGLSFEAACAMPVVAVTMIECFTKARLKAGESILIQTATGGTGLVAVQLAQHAGARIYATAGSAAKLNYLAGLGIEHRINYLEQDFEAELMRLTGGRGVDVVINTLGGDAVQKGLNCLAPEGRYIEIAMTALKSAHAIDLSGLANNQTLHSIDLRKLGRTNPAALERGVREMTRLLEAGVISPVLSRIFDFEQVQDAYRWLEDRRNIGKVVVSVPLTYRYQAPDSGERIAIEPIAVIGMSGRFARAGDLQELWQALAGGEDLIEEASRWPLDALGPDQEPYCHHGGFLRDIDAFDPMFFNISGHEAAVMDPQQRIFLEEAWRALEDAGYAGASVEGRRCGVYVGCAAGDYQRLLERDAPAQAFWGNAGSLIPARIAYHLDLQGPAVAIDTACSSSLVALHQACQALRHGEAELALAGGVFVQSTEHFYLQANRAGMLSPRGRCHTFDARADGFVPGEGAGVVVLKRLSQALADGDLIHGVIKGSGINQDGATNGITAPSARSQTRLEREVYQRHGIDPQQIQVVEAHGTGTVLGDPIEYRALTEALLDGKPTGELGTRCAIGSIKSNLGHTAAAAGIAGVIKLLLALRHRRIPPSLHFEQGNTHIDFSRSPLYVPTTLEDWPASVGGQRLAAVSSFGFSGTNAHAVIGEAPATNRVLPRRPTYLVVLSARTQPQLRLQLERILAHLDGEVPPLASLSHTLLLGRRHFDCRWACLASDLKGLRAQLAEALQKETVGGRIGGADHVPPGEEQLDPLQRRMSDYTDMISHEAARDLLEALREAYLQGLPLDWSFLFQGDGWQRVGLPGYPFARERYWVPVRKTVANTEAALASEPLGQAGMSQGDAGETRTLTWVPQWRSRALVSDAQTAAADRHVVLCELDADPALVEGAVTVRQWTHEGSLDQRYAHYAEYLLTEIQTLAKHRSRRSVLLQLVVPARSEGAVLQALGGLLRTAEKEYSWLRTQLIAVDDVAHLSERLNEEATADPTSRVRYLGSSREVLSYVPAVPACGSAPVWREGGVYWITGGLGGLGLVFAEGIARQVRCPTLVLSGRRAPAPTQQARLERLGELGATVECHALDVSEAVAVAALAQSIVARHGRLSGVIHAAGVLRDGLLHNKRREDLQAVLAPKVAGLLTLERATAGLSLDWLLLCSSVAGVLGNLGQGDYAAANAFLDAYALYRQAPYDRVTRRTRLYSISWPLWEEGGMRIDTDTQAALWREAGVKAMPSEAGLQALNCVLTQDFAHALVMHGDARRLSQVVEGAMPEAVDEEKAATLQSNPADLKAKLEAELAGMIASHLQLPIEALSRDARISEFGFDSISLKAFLKLLNRRHGLALSPAVFFEEPSIRALAAYLLREHGEAFVAIESSAQATLEAPPEPPPVAATSSSTEREQQGGAKPTQREAIAVIGLDGYFPASTDLQEYWDNLWTGRDCITEVPARRWSLDEFYTEDVEVAIRQGRSYSKWGGFIRDIEALDPGFLSGVPAKARQQLNEEQKLFVGIVDRLLVTSGYTEQRLEALRCRRVGVYLGMTAERSAPTDATTSGRNDSPGTLAGMVSRMFRFNGPSVAVDAHSASSMTAVHMACNNLLHSECDAAIAGGVSLLYPDTYRDGCQISLLASHPESRSFSEDKDGVLLADGVGAVLLKRLSTAVEDGDRILAVIRSTVAQSVSSGLSDLPKPELVAASIRENFARATVDPRTISYVEAASAGFPIGDVIEMSATALAFRAYTDQRQFCALGSVKANIGHATAASGISQLAKVVLQLQNGRLAPSIKVGPEQVQAQLRKSPFYVQQQAEDWQRPRLSIDGNEEGREYPRRAMINSMGHGGFYAGAILEEYCGPVLEDQ",
        ],
        "9stFB1fGjCdZVZWHLVI3OD4A_DV3WcV6": [
            "AXA20097.1",
            "AXA20097.1",
            set(),
            "MSAEANKAIVTAMYEALNNRDAKGHFGHMADDVQVTYFGNHRFSRTFHGKQDLFENFTKHFMEYLEGPLDFRVGNIIATDDYAVIEGQGIGRTKDGQDYNNVYCIVMRLVEGKVTEIREYMDTDLAKRIFG",
        ],
        "MSHRSCZfdBJP8vdJdaXfeZrThH_4EUMm": [
            "AXA20098.1",
            "AXA20098.1",
            set(),
            "MKDLTQGAVTRHIVSMAVPIGVGAVFQSLYYLIDLYFVGCLGSDALAGVSAAGNLSLLVMALTQVLGAGTLALMAQAAGRKNEDQARGIFNQALVLSICSGLVLLLLGYALTGVYLRSTSATLVVAEQGQRYLYWFLPGMAFQFVLTAMASALRGIGVVKPTMMVQLMTVVLNIVLAPVLIVGWGTGYAMGVAGAGLASSIALAVAAPAMAWHFHQLGHYVQVRRELLRPRWADWRRLLTIGLPAGGEYTLLFLYSAATYWAIRDFGPAAQAGFGAAARLMQVLLLPALALSFATGPMVGQNLGAGMAERVRRTFGAAIWLISSIMLLACLLLLWEGETLLGRFATGGEVIEQGTHYLHIACWSFIAQGVIFTCSSVFQGLGNTRPALVSSVVRLLCFVLPVAWLSARGHFPIIAVWYLSLASIFLQALLSLWLVRQEFRLKLAPVAGLPSMQAIRR",
        ],
        "Qi23auOUTcBzWTmDHuGinrzIuqH7-zVn": [
            "AXA20099.1",
            "AXA20099.1",
            set(),
            "MSSGYQPDAVQARGNPARGGTVFMFGGQGTQYFQMGRELYRSHPVFRERMDGCDALIRAELGYSLTEVLYEAGHTSSTPFDEILYTHPALFSIGHSLGEAIRADGIEPKAVLGYSLGEYIGLVAAGCLGWEGGLRLVIRQAQVLAQHGAAGGMLSVFAELEQFWQRPDLYTGSQLAVVNFSGSFCVAGAPEAMEAIKRRLDEEQRISTLLPIRFAFHSSGIDALEQRIRGLTADLRIKPGRVPAYSCMLRRAIGTAEFADPTAYCWRVVRDSVHFEQTARRLAAEIPVPFMVDLSATGTLATFIKYARIDGVSYSHCINQFGRDLNEYKSLMQRFAV",
        ],
        "nk3UJUyLaWr8LochtohJ9L6eugdChZL9": [
            "AXA20100.1",
            "AXA20100.1",
            set(),
            "MSSVATVPTTPSYTMRPLLHAVLFLRFRHARVPLIGQRIRLVRRGVMARLMQLPVTVVENRSKTVATTGGSLAENVVYRRFKTLTGNCFWGASHRLAGAQGHRSRRRNQPHGGPRSSIIRSYRLKLCP",
        ],
        "du1Ncfm5UYiFYgDWD8KW1AQJNHlAcVXL": [
            "AXA20101.1",
            "AXA20101.1",
            set(),
            "MSAIPVHASPFLFELLTNALPDVYQSYTELRRLGGVVRDSSGVILLARHDDIVKTFADRRFTSVGRQSLASMSPALQEKMKGSRLFELLNFRDGDSHRQARRVLAGLFNPDTMEWIRVEIAAEAGVLMERWQTGEADDFVAAVTAQLPMRALAMLIGIEEVGLQDFFLRTRNFGAWLSASVFSEEGLSLVADEFVWARGWLRAQLVGKPLFDAVSDDTRENLLGDLLLVLVTGYDSSVALLGNGLATLVSVPALRDRLARDAGLSLRAAEELIRYDSPAQVAFRLALEDVTLGEHRIAKGEFVALVNGSGNRDETIYEHADQLNFDRPKQRALAFGSGPHACAGAALAKAQLTGFLDGVRPWLPHLTLDDATAPPSQHGLLRYRTHLMLRHTVV",
        ],
        "5IYMhENey2WCMrKPUz3AqBIZuFSv6DPP": [
            "AXA20102.1",
            "AXA20102.1",
            set(),
            "MAYRSTIYLDWGTFFLIETPAASIQPLLPASVRPLLAASGAAILTLNVVHFLAGGEQVDLPANHEIDIGVLVELDNSEFEADLPQASVAIHVLKVASTARGYVDLCQNTGYRVIDPVALAFEINPATLTASVCDADGPILGCRQLDLDLEYDEFRRIGQDVMYDAERGIHRANYIFSGKGLSRPMTNAMELRVHAHPFFADLGIEPGVLVCLDQFALKPSSRSSLAFYRPNQVEMDETVQAEKD",
        ],
        "iiWqYfcbDGjauCrUsdiI1pAlG5Syx_-L": [
            "AXA20103.1",
            "AXA20103.1",
            set(),
            "MSTRMFITGGACGLGRALAQRYARDGASVCIGDVDDAAGRQVVEALLAEGATAHYLHCDVTREPQLHEAANWLRTQWGGVDIVINNAGVAQIGGITESSLEDWQWAVDINLLGVVRGCKAFIPLLQAQGGGKLLNVASIAGLLYMPKSGGYNATKAAVIALSETLQLELHDSGIQVSVACPAYFRTDLARNMRASNAQLRQRTHNLVEQARLGSQEVAELIHAGLARGDTHILTHPATRIAWRLKRWLPYRWYLGIARKQIAKAGTAAEDAPA",
        ],
        "WbViYzQw8y-XfCQMgQXkedGduNMJPa14": [
            "AXA20104.1",
            "AXA20104.1",
            set(),
            "MSFTLQGIPVSRGISIGRAYLIAPAALDVAHYLVEAQQLEAEIERLRAALLAVRGELDVLRSDLTEDTPSEVAAFIDVHSMILNDALLVQETIDLIRVRRYNAEWALTKQLDVISGHFDDIEDEYLRERKADIEQVVERVLKALAGEPSASQALGGATGGKDEMIVVAHDIAPTDMIQFKNQAFHAFVTNLGGHTSHTAIVARSLGIPALVGVQHASALIRQNDLIIVDGERGIVIVDPAPIVLEEYSYRQSEMALEQHKLQRLKFSPTQTLCGTQIDLMANIELPDDTKTAVEAGAVGVGLFRTEFLFMNDANALAEEEEQFRAYRRTVEMMNGKSVTIRTIDVGADKPREEGYETAPNPALGLRAIRWSLSAPQMFLTQVRAILRASSVGPVKILVPMLAHAQEIDQTLNLIREAKRQLDDAGLAYDPNVRVGAMIEIPAAAIALPLFLKRVDFLSIGTNDLIQYTLAIDRADNTVAHLYDPLHPAVLHLISGTLRAAKRAGVPVSVCGEMAGDPAMTRLLLGMGLTEFSMHPSQLLVVKQEILRAHLKAIEKSTADVLASYEPEDVQLALKQLAAAKPKADLAA",
        ],
        "4_8182J88axMDpFJBZI6kLNJAu8Ittm3": [
            "AXA20105.1",
            "AXA20105.1",
            set(),
            "MIVPMKNSDLSLPSKLDAAILRGRDALAQRQNADGSWRFELESDATITAEYILMMHFIGKIDDVRQARMARYLREIQRLAAHGGWDLYVDGAPDVSASVKAYFALKAAGDSEDAPHMARARETILHLGGAARCNVFTRILLVTFGQVPWRATPFMPIEFVLFPKWVPISIYKVAYWARTTMVPLLVLCSLKARAKNPRGISIRELFVTAPEEERNYFARGGFIRNLFLYIDRTLCSLDALIPKALRRRAIRHAESWCAERMNGEDGMGGIFPPIVYSYQMMEVLGYAEDHPLRRACEDALEKLLVERPDGSMYCQPCLSPVWDTAWATMALEQARAVPDTREQPTVSAAQLQVGITRACDWLAGKQVTELKGDWIENAPTETPAGGWAFQYENPYYPDVDDSAVVAAMLHQRGCAMARLTGTDPYTEVVSRGLDWIRGLQSRNGGFGAFDADCDRLYLNLIPFADHGALLDPPTEDVSGRVLLCFGVTGRAEDRSALARAIEYVKRTQREDGSWWGRWGTNYIYGTWSVLAGLALAGEHCSQPYIARAIDWLCARQNADGGWGETNDSYVDPSLAGTNGGESASNFTAWALLAQMAFGEWQSESVQRGIRYLLSVQQADGFWWHRSHNAPGFPRIYYLKYHGYTAYFPLWALARYRCLSQAHVATSASPAAEARRGVVL",
        ],
        "IRqRpDzrGB9UhHJD6AzDq_6Xupj00Nte": [
            "AXA20106.1",
            "AXA20106.1",
            set(),
            "MNFDDYCQKKAAPPGSSIYYALRQAPLTSQGALIALFALRRELEEATRETSEPAIGQTKLAWWRKELAALAEGQPSHPVSQALALHASAIASDHAVLQALADGFAMDLEQTRYLDWPNLRRYVERVGGGFAGLVARATTGPAIPADIAPAWAASLGEGLQLAQIVEDVGDDLRHGRVYLPFDELQRYEVTSVDLMHRRYSPAFRKLMRFQVTRARETLHAALEAIPAAELPAQRSLRALAALALATLDEIEGEDYQVLHQRILLTPIRKLWVAWRAAQRRY",
        ],
        "ewjVum5PbpEJA4rl-BfnCAypKl5HXb7x": [
            "AXA20107.1",
            "AXA20107.1",
            set(),
            "MASMGIISGNTNGLRVSEVRHTVDDADANGDLGRLRGDVTRPETVAGEGLEPIHRILGKRSPVVATFLLPFSTTVTGNCINRAIMPRRTGHIRWPMNDTLAWRNRRNSTACSNGRMAWLGVVGTITANDIDLFFAWNLVEQLGQSITVSHILIRH",
        ],
    }
    assert {k: v.sequence for k, v in sg_object.proteins.items()} == PROTEIN_DICT


def test_fasta_string_parse():
    sg_object = SocialGene()
    sg_object.parse(">asdads\n dasfa")
    assert list(sg_object.proteins.keys())[0] == "lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ"
    protein_parse_results = [
        {k: [v.description, v.external_id, v.domains]}
        for k, v in sg_object.proteins.items()
    ]
    assert protein_parse_results == [
        {"lZA5w9NxutVBhoqjqfsX_GX_dfugQZLQ": ["asdads", "asdads", set()]}
    ]
