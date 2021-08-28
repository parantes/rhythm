# Barbosa's dynamical rhythm model
#
# author: Pablo Arantes <pabloarantes@protonmail.com>
# created: 2019-01-18
#
# R implementation of Barbosa's coupled-oscillator model of
# rhythm production. Original Matlab source code can be found in
# Barbosa (2006, pp. 480-485).
#
# Variables
# ---------
# alpha: entrainment rate
# beta: decay rate
# T0: syllable oscillator resting period
# w0: relative coupling strength
# catalexis: remaining syllables in the last stress group
# resetting: syllable-oscillator period resetting method ("fixed" or "variable")

rhythm <- function(
  alpha = 0.4,
  beta = 1.1,
  w0 = 0.78,
  T0 = 0.165,
  units_per_sg = c(4, 4),
  amplitudes = c(0.5, 1),
  catalexis = 0,
  resetting = "fixed",
  plot.contour = F
  ) {
  
  # Check if number of stress groups and amplitudes are the same
  if (length(units_per_sg) != length(amplitudes)) {
    stop("Number of stress groups and amplitudes must be the same.")
  }

  # Vector to hold synchronicity function values
  sync <- numeric(sum(units_per_sg))
  # Number of stress groups
  nsg <- length(units_per_sg)
    
  # Compute the synchronicity function for every VV unit
  # in the simulated utterance
  # ----------------------------------------------------
  
  # Count the number of VV units along utterance
  vv_count <- 1

    for (sg in 1:nsg) {
    vvs <- units_per_sg[sg]
    for (vv in 1:vvs) {
      # Sync function
      if (vv == 1) {
        # 1st unit
        s <-  w0 * exp(-vvs + 2)
      } else if (vv == vvs) {
        # Last unit
        s <-  w0 * exp(-5.81 + (0.016 * T0 * 1000))
      } else {
        # General case
        ## units are counted beginning at 0, not 1 -
        ## that's the reason to 'vv-1'
        s <-  (1 - w0) * sync[vv - 1] + w0 * exp(-vvs + (vv - 1) + 2)
      }
        # Append computed 's' value to the 'sync' vector
        sync[vv_count] <- s
        # Update counter
        vv_count <- vv_count + 1 
    }
  }
  
  # Sync value for VV units in catalexis are equal to
  # the last computed sync value
  if (catalexis > 0) {
    for (cat in 1:catalexis){
      sync <- append(sync, tail(sync, 1))
    }
  }
  
  # Compute entrained syllable oscillator abstract period duration
  # --------------------------------------------------------------

  # Vector to hold simulated duration values
  sim_dur <- numeric(length(sync))
 
  if (catalexis > 0) {
    # Append catalexis units as a stress group
    units_per_sg <- append(units_per_sg, catalexis)
    # Amplitude of phrasal stress in catalexis is a copy
    # of the last one
    amplitudes <- append(amplitudes, tail(amplitudes, 1))
  } 
  
  # Count the number of VV units along utterance
  vv_count <- 1

  for (sg in 1:(nsg + catalexis)) {
  # Number of syllable-oscillator periods over which
  # period resetting will act 
    if (resetting == "fixed") {
      res_len <- 2
    } else if (resetting == "variable") {
      res_len <- round(0.7 * units_per_sg[sg])
    } else {
      stop("Resetting method has to be \'fixed\' or \'variable\'.")
    }
    # Number of VV units in the current SG
    vvs <- units_per_sg[sg]
    for (vv in 1:vvs) {
      if (vv_count == 1) {
        # Utterance-initial unit is special
        if (w0 == 0){
          dur_prev <- T0
        } else {
          dur_prev <- 7 * (T0^2)
        }
      } else {
        dur_prev <- sim_dur[vv_count - 1]
      }
      s <- sync[vv_count]
      amp <- amplitudes[sg]
      # Period entrainment computation
      deltaT = (alpha * dur_prev * s * amp)
      # Reset component
      if (sg == 1){
        # For the SG in initial position, period resetting
        # is not modulated by previous phrase stress amplitude
        amp_prev <- 1
      } else {
        amp_prev <- amplitudes[sg - 1]
      }
      reset <- -beta * (dur_prev - T0) * amp_prev
      # Reset is applied only to VV units periods that are
      # are not in catalexis and are within the reset window  
      if ((sg <= nsg) & (vv <= res_len)) {
        deltaT <- deltaT + reset
      }
      dur <- dur_prev + deltaT
      sim_dur[vv_count] <- dur
      vv_count <- vv_count + 1
    }
  }
  if (plot.contour == T) {
    plot(cumsum(sim_dur),
         sim_dur,
         type = "b",
         xlab = "time (s)",
         ylab = "abstract duration (s)"
    )
  }
  return(sim_dur)
}
