function [dprime, d_primes, all_possible_rates] = calc_dprime(hit_rate, false_alarm_rate)
% Calculate d-prime based on hit rate and false alarm rate

% Calculate d-prime (% use 0.001 and 0.999 as z-values cannot be calculated
% for 0 and 1 (-inf and inf, respectively))
all_possible_rates = [0.001 0.01:0.01:0.99 0.999];
[all_false_alarms, all_hits] = meshgrid(all_possible_rates, ...
    all_possible_rates);
zvalue_false_alarms = norminv(all_false_alarms);
zvalue_hits = norminv(all_hits);
d_primes = zvalue_hits - zvalue_false_alarms;
dprime = d_primes(dsearchn(all_possible_rates', hit_rate), ...
    dsearchn(all_possible_rates', false_alarm_rate));
% a_primes = sign(all_hits-all_false_alarms) .* (((all_hits-all_false_alarms).^2 ...
%     + abs(all_hits-all_false_alarms)) ./ (4*max(all_hits,all_false_alarms) ...
%     - 4 * all_hits.*all_false_alarms));
% aprime = a_primes(dsearchn(all_possible_rates', hit_rate), ...
%     dsearchn(all_possible_rates', false_alarm_rate));