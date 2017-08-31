p = Array.new
while line=gets
  next if line=~/^\#/
  a = line.split(/\s/)
  t = a[0].to_f
  next if t < 20.0
  p.push a[2].to_f
end

puts p.inject(0.0){|r,i| r+=i}/p.size
