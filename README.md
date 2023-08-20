# Evolutionary_RPS
Implementation of three species May-Leonard model with evolution


Code Properties
1. Evolution is performed by sampling a normal distribution about the parent's predation probability. If the number chosen will result in an invalid probability, a new number is chosen.
2. The other actions have their probabilities changed based off of what happens with predation (Ex. If predation increased by 0.2, diffusion, birth, swapping, and death will decrease by 0.05) in order to preserve a total probability of 1
3. The particles live on "sheets". The code forces all particles to the lowest sheet possible meaning that there can be no gaps in between particles between sheets of the same x,y coordinates. This was done in hopes of reducing search time during predation and swapping. However, this could induce a bias towards particles on the lowest level.
