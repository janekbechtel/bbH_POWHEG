# two stages of importance sampling grid calculation
echo "start"
cat powheg.inputSave > powheg.input
echo "parallelstage 1" >> powheg.input
echo "xgriditeration 1" >> powheg.input
cp powheg.input powheg.input-1-1          # save a copy of the
# powheg.input for this stage
for i in {1..12}
do
echo $i | ../../pwhg_main > run-st1-xg1-$i.log 2>&1 &
done
wait
cat powheg.input-save > powheg.input
echo "parallelstage 1" >> powheg.input
echo "xgriditeration 2" >> powheg.input
cp powheg.input powheg.input-1-2          # save a copy of the
# powheg.input for this stage
for i in {1..12}
do
echo $i | ../../pwhg_main > run-st1-xg2-$i.log 2>&1 &
done
wait
# compute NLO and upper bounding envelope
# for the generation of the underlying born configurations
cat powheg.input-save > powheg.input
echo "parallelstage 2" >> powheg.input
cp powheg.input powheg.input-2          # save a copy of the
# powheg.input for this stage
for i in {1..12}
do
echo $i | ../../pwhg_main > run-st2-$i.log 2>&1 &
done
wait
# compute upper bounding coefficients for radiation
cat powheg.input-save > powheg.input
echo "parallelstage 3" >> powheg.input
cp powheg.input powheg.input-3          # save a copy of the
# powheg.input for this stage
for i in {1..12}
do
echo $i | ../../pwhg_main > run-st3-$i.log 2>&1 &
done
wait
# generate events
cat powheg.input-save > powheg.input
echo "parallelstage 4" >> powheg.input
cp powheg.input powheg.input-4          # save a copy of the
# powheg.input for this stage
for i in {1..12}
do
echo $i | ../../pwhg_main > run-st4-$i.log 2>&1 &
done
wait
