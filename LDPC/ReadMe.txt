Connessione rete unitn:
globalprotect connect --portal vpn.icts.unitn.it --username giovanni.tognolini@unitn.it

Connessione server unitn:
ssh tognolini@10.194.32.83
psw: ----

Connessione server Ancona
ssh root@193.205.130.203
psw: pippopippo


scp -r /home/triki/Scrivania/multiKey/ tognolini@10.194.32.83:/home/tognolini/
scp -r /home/triki/Documenti root@localhost:/home/Sagemath