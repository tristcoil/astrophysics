# Sobral 1919: station, weather, timing, comparison plates and plate selection

Source: Dyson, F. W., Eddington, A. S. & Davidson, C. 1920, *Phil. Trans. R. Soc. A* 220, 291-333, "A Determination
of the Deflection of Light by the Sun's Gravitational Field, from Observations made at the Total Eclipse of May 29,
1919". Section III, "The Expedition to Sobral" (observers Dr. A. C. D. Crommelin and Mr. C. Davidson), pp.296-302,
plus the passages on pp.306, 309-312 and 330-331 that say which plates were used and why.

Scans: archive.org item `philtrans06337895`, leaves 4-11, 15, 18-21, 39-40 (leaf N = printed page 291 + N), read
on the full-resolution JPEG 2000 files. Quotations keep the printed wording and spelling ("coelostat" is printed with
the ligature). The raised decimal point is written as a full stop (74°·4 becomes 74°.4, 1″·98 becomes 1″.98), and
"9/10" is printed as a stacked fraction. Times are **astronomical G.M.T.**, counted from noon, as printed. So
"May 28, 23h. 58m. 18s." is civil 1919 May 29, 11:58:18 G.M.T. Values marked *derived* are not in the paper.

The exposure tables of p.300 are in `dyson1920_sobral_exposures.csv`. Checks and discrepancies are in
`dyson1920_sobral4inch_results_notes.md`.

## 1. Station (p.296, §9; p.297, §10)

> "9. Sobral is the second town of the State of Ceara, in the north of Brazil. Its geographical co-ordinates are:
> longitude 2h. 47m. 25s. west; latitude 3° 41′ 33″ south; altitude 230 feet. Its climate is dry and though hot not
> unhealthy." (p.296)

> "A convenient site for the eclipse station offered itself just in front of the house; this was the race-course of
> the Jockey Club, ..." "We laid down a meridian line, after which brick piers were constructed for the coelostats
> and for the steel tube of the astrographic telescope." (p.297)

**The printed longitude is almost certainly wrong** (notes, R6). 2h 47m 25s W is 41.854° W, about 167 km west of
Sobral. The project uses 40.35° W, 3.69° S for the station (`fetch_data.py`), which is 2h 41m 24s W. That is one digit
away: "41m" printed as "47m". The paper's own timing confirms this (*derived*, JPL Horizons, topocentric Moon and Sun
at 5 s steps, latitude as printed):

| Station longitude | Second contact (UT) | Third contact (UT) | Totality |
|---|---|---|---|
| 40.35° W (2h 41m 24s, the race-course) | 11:58:10 | 12:03:28 | 5 m 18 s |
| 41.854° W (2h 47m 25s, as printed; run at 41.857° W, see Audit) | 11:55:54 | 12:00:21 | 4 m 27 s |

The paper assumes second contact at 11:58:18. Its exposures run from 11:58:21 to 12:03:23, and the 4-inch plates
exposed after 12:00:21 show stars. Only the 40.35° W station fits.

The field altitudes printed for the July comparison plates (p.301) imply 2h 42m 48s to 2h 44m 06s W (*derived*).
That is closer to 41m than to 47m, with about 1 minute of uncertainty from the field centre and from which part of
the exposure the altitude refers to. The printed latitude, 3° 41′ 33″ S (3.6925° S), is 0.3 km from the project's
3.69° S. 230 feet is 70 m.

## 2. Instruments and focusing (pp.297-298, §10)

- **Astrographic** (13-inch object glass of the Greenwich astrographic, in a steel tube, fed by a 16-inch
  coelostat):
  > "for convenience it was adjusted so that the directions of R.A. and declination were parallel to the sides of the
  > plate; this involved a tilt of the plate holders of about 4 degrees to the horizontal." (p.297)

  > "The focussing was at first done visually on Arcturus, using an eyepiece fitted with a cobalt glass (after the
  > plate supports and object-glass had been adjusted for perpendicularity to the axis). A series of exposures was then
  > made, the focus being varied slightly so as to cover a sufficient range. Examination of these photographs showed at
  > once that there was serious astigmatism due to the figure of the mirror of the 16-inch coelostat. By inserting an
  > 8-inch stop this was reduced to a large extent, and this stop was henceforth used throughout; but the defect was of
  > such a character that it was clear that it would be necessary to stay at Sobral and obtain comparison plates of the
  > eclipse field in July when the sun had moved away." (p.298)
- **4-inch** (p.295: Father Cortie suggested "the use of the 4-inch telescope of 19-feet focus, which he had used at
  Hernosand, Sweden, in 1914, in conjunction with an 8-inch coelostat, the property of the Royal Irish Academy"; p.296:
  "With the 4-inch lens Father CORTIE lent the square wooden tube used by him in 1914"; the 28 s exposures are
  mentioned on p.297):
  > "The 4-inch lens was taken as an auxiliary; we used the square wooden tube, 19 feet in length, originally used by
  > Father CORTIE at Hernosand in 1914, together with the 10 × 8-inch plate carriers. Study of the star-diagram showed
  > that seven stars could be photographed by turning the plate through 45 degrees." (pp.297-298)

  > "The 16-inch coelostat had free adjustment for all latitudes; but the 8-inch one, constructed for European
  > latitudes, was mounted on a wooden base, inclined at an angle of about 40 degrees, constructed before leaving
  > Greenwich. The clock had to be separated from the coelostat, mounted on a wooden base and reversed, to adjust to
  > the Southern Hemisphere. It performed very satisfactorily, and no elongation of the star images is shown with 28
  > seconds' exposure." (p.297)

  > "The focus of the 4-inch was determined in a similar manner. The images, though superior to those of the
  > astrographic, were not quite perfect, and here again comparison plates in July were necessary. Once the focus had
  > been decided on, the breech end was securely screwed up to avoid any chance of subsequent movement." (p.298)
- > "A few check plates of the field near Arcturus were taken, but have not been used." (p.298)
- Driving: "The performance of the 16-inch coelostat was unsatisfactory in respect of driving. There was a clearly
  marked oscillation of the images on the screen in a period of about 30 seconds. For this reason exposure time was
  shortened, so as to multiply the number of exposures in the hope that some would be near the stationary points."
  (p.299)

## 3. Weather during the stay, May 1919 (p.298, §11)

> "The barometer record was interesting in that it showed very little change from day to day, in spite of changes in
> the type of weather; there was, however, a very well marked semi-diurnal variation, with range of about 0.15 inch.
> The temperature range was fairly uniform, from a maximum of about 97° F. towards 3 p.m. to a minimum of about 75° F.
> at 5 a.m. The relative humidity (as shown by a hygrograph belonging to the Brazilian Commission) followed the
> temperature closely, varying from 30 per cent. in the afternoon to 90 per cent. in the early morning."

> "May is normally the last month of the rainy season at Sobral, but this year the rainfall was very scanty; there
> were a few afternoon showers, each ushered in by a violent gust of wind; and on May 25 there was very heavy rain, ...
> There was a fair amount of cloud in the mornings, but the afternoons and nights were clear in the majority of
> cases."

No absolute barometer reading, and no temperature or pressure used for the refraction corrections of p.307, is
printed. Plates were developed at night or before dawn, with water cooled "to about 75°", and "Formalin was used in
every case to harden the films" (pp.298-299). One of the two plate brands proved "unsuitable in the hot climate"
(p.299).

## 4. The eclipse, 1919 May 29 (pp.299-300, §§13-14)

> "13. The morning of the eclipse day was rather more cloudy than the average, and the proportion of cloud was
> estimated at 9/10 at the time of first contact, when the sun was invisible; it appeared a few seconds later showing a
> very small encroachment of the moon, and there were various short intervals of sunshine during the partial phase ...
> As totality approached, the proportion of cloud diminished, and a large clear space reached the sun about one minute
> before second contact. Warnings were given 58s., 22s. and 12s. before second contact by observing the length of the
> disappearing crescent on the ground glass. When the crescent disappeared the word "go" was called and a metronome
> was started by Dr. LEOCADIO, who called out every tenth beat during totality, and the exposure times were recorded
> in terms of these beats. It beat 320 times in 310 seconds; allowance has been made for this rate in the recorded
> times. The programme arranged was carried out successfully, 19 plates being exposed in the astrographic telescope
> with alternate exposures of 5 and 10 seconds, and eight in the 4-inch camera with a uniform exposure of 28 seconds.
> The region round the sun was free from cloud, except for an interval of about a minute near the middle of totality
> when it was veiled by thin cloud, which prevented the photography of stars, though the inner corona remained visible
> to the eye and the plates exposed at this time show it and the large prominence excellently defined. The plates
> remained in their holders until development, which was carried out in convenient batches during the night hours of
> the following days, being completed by June 5." (p.299)

> "14. No observation of contact times was made, but it is known that these times were somewhat before those
> calculated. As the times recorded were reckoned from second contact, it is assumed that this occurred May 28, 23h.
> 58m. 18s. G.M.T." (pp.299-300)

- Wind: "Happily dead calm prevailed at the critical time." (p.299)
- Exposures (p.300, in the CSV): astrographic, 19 plates from 28d 23h 58m 23s to 29d 0h 3m 18s (5 s and 10 s
  alternately). 4-inch, 8 plates of 28 s from 28d 23h 58m 21s to 29d 0h 2m 52s. The last exposure ends 302 s (4-inch)
  and 305 s (astrographic) after the assumed second contact. All exposures fall inside the 5 m 18 s totality computed
  for the station (§1, *derived*).
- > "With the astrographic telescope 12 stars are shown on a number of plates, and seven stars on all but three (Nos.
  > 13, 14 and 19). Of the eight plates taken with the 4-inch lens, seven show seven stars, but No. 6, which was taken
  > through cloud, does not show any." (p.300)

  The plates with no stars or fewer stars date the cloudy minute (*derived*): 4-inch No. 6 (29d 0h 1m 34s to
  0h 2m 2s) and astrographic Nos. 13-16 (0h 1m 36s to 0h 2m 35s). That is 3¼ to 4¼ minutes after second contact,
  a little after mid-totality. p.310 names Nos. 15 and 16 as "taken through clouds" (§8).

## 5. Temperature during the eclipse, from Dr. Morize (pp.300-301)

> "The following table of temperatures, communicated by Dr. MORIZE, and converted into the Fahrenheit scale, shows how
> slight the fall was during totality, probably owing to the large amount of cloud in the earlier stages which checked
> the usual daily rise." (p.300)

The table on p.301 is printed in four blocks of three rows, headed "G.M.T." (d. h. m.) and "Ther." (°). Here it is as
one list, with day and hour carried down where the page leaves them blank:

| G.M.T. (astronomical) d h m | Ther. °F | *derived*: °C | *derived*: civil G.M.T. May 29 |
|---|---|---|---|
| 28 22 45 | 82.4 | 28.0 | 10:45 |
| 28 23 0 | 84.2 | 29.0 | 11:00 |
| 28 23 15 | 82.4 | 28.0 | 11:15 |
| 28 23 30 | 80.6 | 27.0 | 11:30 |
| 28 23 45 | 82.4 | 28.0 | 11:45 |
| 29 0 0 | 80.6 | 27.0 | 12:00 |
| 29 0 15 | 82.0 | 27.8 | 12:15 |
| 29 0 30 | 82.4 | 28.0 | 12:30 |
| 29 0 45 | 83.1 | 28.4 | 12:45 |
| 29 1 0 | 83.8 | 28.8 | 13:00 |
| 29 1 15 | 84.2 | 29.0 | 13:15 |
| 29 1 30 | 84.2 | 29.0 | 13:30 |

Every °F value is a whole or tenth-degree Celsius reading converted, which is consistent with "converted into the
Fahrenheit scale". Totality (11:58-12:03 civil) falls at the 29d 0h 0m reading, 80.6 °F. That equals the 23h 30m
reading and is 3.6 °F below the 23h 0m reading.

## 6. The comparison plates, July 1919 (pp.301-302, §15)

> "15. On June 7, having completed the development, we left Sobral for Fortaleza, returning on July 9 for the purpose
> of securing comparison plates of the eclipse field."

> "Before our departure we dismounted the mirrors and driving clocks which were brought into the house to avoid the
> exposure to dust. The telescopes and coelostats were left in situ. Before removing the mirrors we marked their
> positions in their cells so that they could be replaced in exactly the same position."

> "After our return to Sobral the mirrors and clocks were remounted; the photography of the eclipse field was
> commenced on the morning of July 11 (civil). The difficulty of finding the field with the coelostats was overcome by
> making a rough hour circle on the heads of the coelostats out of millimetre paper." (p.301)

The p.301 table "of exposures made on the field for comparison with the eclipse photographs" is transcribed here as
printed. Date and G.M.T. are astronomical; the reference numbers "follow the civil dates". "Duration." is in s and
"Altitude." in degrees.

| Telescope | Ref. No. | Date | G.M.T. h m | No. of exposures | Duration s | Altitude ° |
|---|---|---|---|---|---|---|
| Astrographic | 11_1 | July 10 | 20 5 | 3 | 5 | 28.9 |
| Astrographic | 11_2 | | 20 16 | 2 | 5 | 31.1 |
| Astrographic | 11_3 | | 20 21 | 1 | 5 | 32.2 |
| Astrographic | 14_1 | July 13 | 20 13 | 3 | 5 | 33.7 |
| Astrographic | 14_2 | | 20 17 | 2 | 5 | 34.5 |
| Astrographic | 14_3 | | 20 19 | 2 | 5 | 34.9 |
| Astrographic | 15_1 | July 14 | 20 15 | 3 | 5 | 34.9 |
| Astrographic | 15_2 | | 20 20 | 2 | 5 | 36.1 |
| Astrographic | 15_3 | | 20 23 | 2 | 5 | 36.6 |
| Astrographic | 17_1 | July 16 | 20 2 | 4 | 3 | 33.8 |
| Astrographic | 17_2 | | 20 15 | 3 | 3 | 36.6 |
| Astrographic | 17_3 | | 20 23 | 2 | 3 | 38.3 |
| Astrographic | 17_4 | | 20 25 | 2 | 5 | 38.8 |
| Astrographic | 18_1 | July 17 | 19 50 | 3 | 4 | 32.8 |
| Astrographic | 18_2 | | 20 1 | 2 | 4 | 34.4 |
| Astrographic | 18_3 | | 20 20 | 3 | 4 | 38.6 |
| Astrographic | 18_4 | | 20 25 | 2 | 3 | 39.5 |
| 4-inch | 14_1 | July 13 | 20 7 | 2 | 25 | 32.4 |
| 4-inch | 14_2 | | 20 16 | 2 | 20 | 34.3 |
| 4-inch | 15_1 | July 14 | 20 17 | 2 | 20 | 35.4 |
| 4-inch | 15_2 | | 20 22 | 2 | 20 | 36.4 |
| 4-inch | 17_1 | July 16 | 20 6 | 3 | 15 | 34.7 |
| 4-inch | 17_2 | | 20 24 | 2 | 15 | 38.6 |
| 4-inch | 18_1* | July 17 | 19 57 | 3 | 20 | 33.6 |
| 4-inch | 18_2 | | 20 24 | 2 | 20 | 39.2 |

The table footnotes read: "The reference numbers follow the civil dates." and "* The 4-inch plate, No. 18_1, was
taken through the glass (see § 17, infra) to facilitate the measurement, and is referred to as the scale plate."

*Derived*: 20h astronomical on July 13 is about 08h civil G.M.T. on July 14, about 05h 19m local mean time at
Sobral (40.35° W; 05h 13m at the printed longitude), before dawn.

> "Thermometer readings, July 10, 74°.4; July 13, 73°.7; July 14, 71°.9; July 16, 72°.3; July 17, 72°.3." (p.302)

These are the same astronomical dates as the plates; the time of day of the readings is not stated. *Derived*: they
are 6.2-8.7 °F below the 80.6 °F read at totality (§5).

> "By July 18 we had obtained a sufficient number of reference photographs. Dismantling of the instruments was
> commenced, and the packing was completed on July 21." (p.302)

The refraction and aberration coefficients of p.307 for the scale plate and the six 4-inch comparison plates
decrease steadily as the altitude printed here increases (notes, check 4d).

## 7. Which 4-inch plates were used, and why (p.302, §§16-17; p.306)

> "16. These photographs were taken on 10 × 8-inch plates. By suitably mounting the camera it was made possible to
> obtain seven stars on the photographs, viz., Nos. 2, 3, 4, 5, 6, 10 and 11 of the table in § 5. Of the eight
> photographs taken during the eclipse seven gave measurable images of these stars, the other plate (No. 6) taken
> through cloud only showing a picture of the prominences."

> "Plates of the same field taken under nearly similar conditions as regards altitude were taken on July 14, 15, 17
> and 18 (civil date). Of these photographs, the second taken on July 14 with two exposures (referred to as 14_2a and
> 14_2b), two photographs taken on July 15 (referred to as 15_1 and 15_2), two on July 17 (17_1 and 17_2), and the
> second photograph on July 18 (18_2) were measured for comparison with the eclipse plates."

> "17. The micrometer at the Royal Observatory is not suitable for the direct comparison of plates of this size. It
> was therefore decided to measure each plate by placing, film to film upon it, another photograph of the same region
> reversed by being taken through the glass. A photograph for this purpose was taken on July 18. This plate is
> regarded merely as an intermediary between the eclipse plates and comparison plates and is referred to as the scale
> plate, being used simply as a scale providing points of reference. ... Each plate was measured in two positions,
> being reversed through 180 degrees, and the accordance of the result showed that the method of measurement was
> entirely satisfactory."

> "The measures, both direct and reversed, were made by two measurers (Mr. DAVIDSON and Mr. FURNER), and the means
> taken. There was no sensible difference between the measurers, which is satisfactory, as it affords evidence of the
> similarity of the images on the eclipse and comparison and scale plates."
> "The value of the micrometer screws (both in R.A. and Decl.) is 6″.25." (p.302)

> "The eclipse and comparison photographs were taken under precisely similar instrumental conditions, but there is the
> difference that the eclipse photographs were taken on the day of May 29, and the comparison photographs on nights
> between July 14 and July 18. A very satisfactory feature of the photographs is the essential similarity of the star
> images on the two sets of photographs." (p.306)

The eclipse plates used are 4-inch Nos. 1-5, 7 and 8, called I-V, VII and VIII in Tables II-VIII. That No. 6 is
"VI" is inferred from the missing VI.

## 8. The astrographic's change of focus, and the decision about it (pp.309-312, 330-331)

> "23. As stated above these photographs were taken with the astrographic object glass stopped down to 8 inches,
> mounted in a steel tube and fed by a 16-inch coelostat. From many years' experience with the object glass at
> Greenwich it is certain that, when the object glass is mounted in a steel tube, the change of scale over a range of
> temperature of 10° F. should be insignificant, and the definition should be very good. It was realised that this
> high standard would not be obtained with the glass used in conjunction with the coelostat taken to Brazil, but
> nevertheless the results shown when the plates were developed were very disappointing. The images were diffused and
> apparently out of focus, although on the night of May 27 the focus was good.* Worse still, this change was
> temporary, for without any change in the adjustments, the instrument had returned to focus when the comparison
> plates were taken in July." (p.309)

> "These changes must be attributed to the effect of the sun's heat on the mirror, but it is difficult to say whether
> this caused a real change of scale in the resulting photographs or merely blurred the images." (p.309)

Footnote, p.309:

> "* The following note made at the time is quoted in full :- "May 30, 3 a.m., four of the astrographic plates were
> developed, and when dry examined. It was found that there had been a serious change of focus, so that, while the
> stars were shown, the definition was spoilt. This change of focus can only be attributed to the unequal expansion of
> the mirror through the sun's heat. The readings of the focussing scale were checked next day, but were found
> unaltered at 11.0 mm. It seems doubtful whether much can be got from these plates.""

> "The photographs were measured in the astrographic duplex micrometer, the eclipse photographs being directly
> compared with the comparison plates taken in July. All the stars shown were measured. They were reduced by the same
> method as that employed for the "4-inch" photographs. With the exception of plates Nos. 15 and 16, taken through
> clouds, the stars numbered 3, 4, 5, 6, 10, 11 and 12 are shown on all the plates; the fainter stars 2, 7, 8 and 9
> are sometimes shown, but No. 1, which is very near the sun, is always drowned in the corona. These plates were only
> measured in declination, as the right ascensions were of little weight." (pp.309-310)

The result and the decision:

- > "Thus the mean value of α obtained from all the astrographic plates is 0″.86, a figure considerably less than
  > that obtained from the 4-inch photographs." (p.310, after Table IX; the screw value there is 1r = 12″.3)
- > "The means of the 16 photographs treated in this manner give α + 243e = + 0r.0435, or with the value of the scale
  > 0r.082 from the previous table α = + 0r.024 = 0″.93 at the limb. It may be noticed that the change of scale
  > arising from differences of refraction and aberration is 0r.020. If this value of e be taken instead of 0r.082 we
  > obtain α = + 0r.039 = + 1″.52 at the sun's limb." (p.312)
- > "The photographs taken with the astrographic telescope support those obtained by the "4-inch" to the extent that
  > they show considerable outward deflection, but for the reasons already given are of much less weight." (p.312)
- > "39. In summarising the results of the two expeditions, the greatest weight must be attached to those obtained
  > with the 4-inch lens at Sobral. From the superiority of the images and the larger scale of the photographs it was
  > recognised that these would prove to be much the most trustworthy. Further, the agreement of the results derived
  > independently from the right ascensions and declinations, and the accordance of the residuals of the individual
  > stars (p. 308) provides a more satisfactory check on the results than was possible for the other instruments."
  > (p.330)
- > "There remain the Sobral astrographic plates which gave the deflection 0″.93 discordant by an amount much beyond
  > the limits of its accidental error. For the reasons already described at length not much weight is attached to
  > this determination." (p.331)

In the paper's words the astrographic result has "much less weight" (p.312) and "not much weight is attached" to it
(p.331). The reasons it gives are the change of focus during the eclipse, attributed to the sun's heat on the
coelostat mirror, and the doubt whether this also changed the scale (p.309). The 0″.93 depends on the scale
difference e = 0r.082 between eclipse and comparison plates being taken from the measures. With the 0r.020 expected
from refraction and aberration alone, the same plates give 1″.52.

The p.300 and p.310 statements about which astrographic plates lack stars differ, though they can be reconciled.
p.300 names Nos. 13, 14 and 19 as the plates without seven stars. p.310 names Nos. 15 and 16, "taken through clouds",
as missing some of the stars 3, 4, 5, 6, 10, 11 and 12. Table IX (p.310) lists eclipse plates 1-12 and 15-18, so
Nos. 13, 14 and 19 were not reduced. It gives No. 15 6 stars and No. 16 7 stars, and p.312 notes that No. 15 does not
show star 10.

## Audit (second reading)

Second, independent reader (T2-audit, 2026-10-20). Every quotation in this file was compared word for word with
full-resolution JPEG 2000 crops of leaves 4-11, 15, 18, 19, 21, 39 and 40. The archive.org OCR text was aligned with
each quotation as a second opinion, and every OCR difference was checked on the scan. All quotations agree with the
print, including every number in them.

- p.296: 2h. 47m. 25s. west (checked at full resolution: "47m"), 3° 41′ 33″ south, 230 feet.
- p.297-299: 40 degrees, 28 seconds' exposure, 4 degrees, 19 feet, 45 degrees, 8-inch stop, 0.15 inch, 97° F at
  3 p.m., 75° F at 5 a.m., 30 and 90 per cent, May 25, about 75°, 30 seconds, 9/10, 58s./22s./12s., 320 beats in
  310 seconds, 19 plates of 5 and 10 s, eight of 28 s, June 5.
- p.300: 23h. 58m. 18s. and the plate key.
- p.302: 6″.25 and the thermometer readings 74°.4, 73°.7, 71°.9, 72°.3, 72°.3.
- p.309-312: 10° F, May 27, 11.0 mm, Nos. 15 and 16, stars 3, 4, 5, 6, 10, 11, 12 and 2, 7, 8, 9, 0″.86,
  1r = 12″.3, 0r.0435, 0r.082, 0r.024, 0″.93, 0r.020, 0r.039, 1″.52.
- p.330-331: §39 and 0″.93.

The p.301 tables were re-read cell by cell with no differences: the 12 temperatures with their times, and the 25 rows
of the comparison-exposure table (date, G.M.T., number of exposures, duration, altitude). The derived values were
recomputed: 3° 41′ 33″ = 3.6925°, 230 ft = 70 m, the °C column, 3.6 °F, 6.2-8.7 °F, and 3¼-4¼ min after second
contact for the cloudy interval.

Corrections made:

1. §1: 2h 47m 25s W is 41.854° W (was 41.857° W). The Horizons run in the table used E-longitude 318.14306°, i.e.
   41.857° W, a 10″ conversion slip. An independent skyfield + DE440s computation gives the same contacts at
   41.854° and 41.857° W to 0.5 s (11:55:56 and 12:00:20.5 UT1), so the tabulated times stand.
2. §1: the station coordinates 40.35° W, 3.69° S are set in `fetch_data.py`; `config.py` does not define them.
3. §2: "Father Cortie's lens" was not the paper's statement. It is replaced by the paper's own words: on p.295, the
   4-inch telescope Father Cortie had used at Hernosand in 1914, with the Royal Irish Academy's 8-inch coelostat; on
   p.296, Cortie lent the square wooden tube.
4. §6: 08h civil G.M.T. is 05h 19m local mean time at 40.35° W (05h 13m at the printed longitude), not 05h 25m.
5. Leaf 4 (p.295) is added to the list of scans.

Independent checks of the *derived* statements (skyfield 1.55 + DE440s, smooth limbs, UT1 with ΔT = 21.2 s,
latitude as printed):

- Totality at 40.35° W runs 11:58:11.5 to 12:03:26.5 UT (5 m 15 s). The table here gives 11:58:10 to 12:03:28
  (5 m 18 s), which agrees to 1.5 s per contact. All exposures fall inside totality.
- At the printed longitude totality runs 11:55:56 to 12:00:20.5 (4 m 24.5 s), against 11:55:54 to 12:00:21 here.
- The p.301 altitudes give zero mean offset at 2h 42m 50s to 2h 44m 04s W. This matches the 2h 42m 48s to 2h 44m 06s
  stated in §1.

Noted, quotation unchanged: p.312 prints "α + 243e" with no decimal point before 243 (checked at full resolution).
The paper's own arithmetic needs 0.243: 0.0435 - 0.243 × 0.082 = 0.0236 → 0r.024, and 0.0435 - 0.243 × 0.020 =
0.0386 → 0r.039.
