sehr schöne Abbildungen! Ein paar Fragen bzw. Kommentare habe ich aber.

Zu S10_new (umgekehrte Invasionsreihenfolge): Ich hatte mir den Fall auch schon grob angesehen und für d_I-Werte, die zu keiner Turing instability führen, einen etwas größeren Bereich gefunden, in denen H_I nicht überlebt, wenn H_S der invader ist. Gegenwärtig habe ich in Fig. S9 eine Zeitserie, die das veranschaulicht, ähnlich Fig. S11 in der vorherigen Manuskriptversion. Wenn deine genaueren Simulationen jetzt zeigen, dass dieser Fall wirklich nur ganz selten, direkt an der Grenze des Koexistenzbereiches auftritt, können wir die lange Erklärung inkl. der Zeitserie auch weg lassen, denke ich.

In S7_new gibt es in beiden panels kritische Bereiche, in denen Koexistenz mit statischer Dynamik angezeigt wird. In panel A ist das der Bereich zwischen den Turing boundaries von H_I, bei sehr hohen a_I. Eigentlich müsste es da Oszillationen geben, da d_S so klein ist, dass H_S eine oTI erzeugt. Könntest du dir zur Verifikation Zeitserien aus dem Bereich anschauen und ggf. das Kriterium zur Feststellung von Oszillationen anpassen? In panel B ist das ganz ähnlich, in dem kleinen Bereich unterhalb der oTI-boundary von H_I. Da müsste es eigentlich auch Oszillationen geben.




Hallo Felix,

jetzt hat es doch sehr viel länger gedauert, als geplant. Ich habe natürlich eine Verlängerung für die Revision beantragt und bin inzwischen auch weit fortgeschritten mit den Überarbeitungen. Ich schicke dir gleich einen Einladungslink zum overleaf-Dokument. Dort habe ich auch die Kommentare der Gutachter und meine Antworten darauf eingepflegt (am Ende des Dokuments).

An dich hätte ich noch folgende Bitten:

1) Die Abbildung, die du mir zuletzt geschickt hast, gefällt Toni und mir sehr gut. Ich habe nur eine Reihe kleiner Änderungswünsche:

i) Bet hedging -> bet-hedging, Maladaptive -> maladaptive
    -> habe ich geändert
    
ii) In den Achsenbeschriftungen von panel A: dS und dI in math font, mit S und I als Index
    -> habe ich geändert

iii) Die Linien für die Turing boundaries müssen noch einmal angepasst werden, weil ich auch in den korrigierten Werten Fehler hatte. Die (hoffentlich jetzt wirklich) korrekten, zweifach unabhängig verifizierten Werte sind (bei a_I = 1.0 und a_S = 1.3): d_I (oTI) = 0.048367, d_I (sTI) = 0.308134, d_S (oTI) = 0.080921, d_S (sTI) = 0.252933.
    -> da ändern sich die boundaries vom inferior competitor

iv) Reviewer 2 hätte gerne Symbole oder ähnliches, die die Zeitserien in Abb. 4 mit dieser Abbildung verknüpfen. Die dispersal rates, die ich in den Zeitserien verwendet habe, sind wie folgt: A: d_S=0.005, d_I=0.018, B : d_S=1, d_I=0.47, C: d_S=0.005, d_I=0.1, D: d_S=1, d_I=0.13. Wenn die weiße Symbole mit schwarzer Umrandung (z.B. Quadrat, Kreis, Raute, Dreieck?) in Abb. 3 einfügst, kann ich entsprechende Symbole in die Zeitserien in Abb. 4 einfügen.
    -> das habe ich gemacht
    
v) Legende: "Coexistence only possible with self-organised pattern formation" bitte um "due to heterogeneity modulation" ergänzen. Außerdem: dynamic -> dynamics
    -> habe ich ergänzt
    -> dynamics habe ich geändert

2) Wir bräuchten eine aktualisierte Version des Flussdiagramms, Abb. S1, mit den vereinfachten nutrient uptake rates und den random (not plastic) heterotroph dispersal rates.
    -> morgen

3) Folgende coexistence plots (analog Abb. 3A) bräuchten wir für die Sensitivitätsanalyse:

i) Standardmodel mit H_I als resident und H_S als invader (analog Abb. S10A im alten Manuskript)
    -> simulation neu erstellen

ii) Koexistenz für H_I mit plastic dispersal (k_I = 2), analog Abb. 2B des alten Manuskripts. H_S bleibt bei random dispersal. Die zu variierenden Parameter sind entsprechend d_S auf der x-Achse (gleicher Bereich wie in Abb. 3) und d_max,I auf der y-Achse (Werte gegenüber Abb. 3 verdoppelt, Position der Turing boundaries entsprechend anpassen!). Für diese Abbildung wäre es ganz gut, wenn du auch die Version mit environmental heterogeneity zur Gegenüberstellung rechnen könntest. Diese Abbildung wäre das Pendant zu Abb. S6 im alten Manuskript, würde aber nur noch aus einem panel bestehen (bzw. einem zweiten, kleinen für die Koex. unter environmental heterogeneity). Die Variante mit plastic dispersal von H_S müssen wir nicht explizit zeigen sondern können sie auch rein argumentativ abhandeln.
    -> morgen programmieren
    
iii) Analog zu Abb. S7 des alten Manuskripts wären zwei Koexistenz-Abbildungen mit a_I (von 0 bis 1.3) auf der x-Achse und d_I (Werte wie in Abb. 3) auf der y-Achse, einmal mit d_S = 0.005 und einmal mit d_S = 1 gut.

Kannst du mir diese Abbildungen im Laufe der nächsten Woche schicken?

Viele Grüße,

Christian