#!/usr/bin/ruby

require 'set'

def calc_edges_count(graph)
    type, size = graph.split('_', 2)
    case type
    when 'rmat'
        pow, degree = size.split('_').map(&:to_i)
        return (2 ** pow) * degree / 2
    when 'ssca2'
        size = size.to_i
        case size
        when 22
            return 134_637_804
        end
    else
    end
    return 0
end

def parse_times_from_file(filename)
    # get all times values
    times = File.open(filename, 'r').each_line.select{ |s| s =~ /^\d+\.\d{3}\s*$/ }.map{ |s| s.chomp.to_f }
    # skip reference time
    times = times[2..-1]
    # convert to pairs
    times.each_slice(2).map{ |pt| { :pre => pt[0], :run => pt[1] } }
end

def best_times_from_directory(dirname)
    Dir.foreach(dirname).reduce(nil) do |memo, runname|
        # skip '.' and '..'
        next memo if runname.start_with?('.')
        # get times for this run
        times = parse_times_from_file(dirname + '/' + runname + '/errors')
        # exclude prepare time
        times = times.map{ |t| t[:run] }
        # set min values
        unless memo.nil?
            times = memo.zip(times).map{ |a| a.min }.dup
        end
        times
    end
end

def generate_data(implementation, graph)
    name = "#{implementation}-#{graph}"
    times = best_times_from_directory(name + '.raw')
    # print times in sec
    File.open(name + '.txt', 'w') do |f|
        times.each{ |t| f.puts t }
    end

    # print times in pgfplot format
    File.open(name + '_time_pgf.txt', 'w') do |f|
        5.times do |pw|
            tcount = 1 << pw
            time = times[pw]
            f.print "(#{tcount}, #{time}) "
        end
    end

    # print performance in mteps
    edges_count = calc_edges_count(graph)
    puts "edges count for graph #{graph} is #{edges_count}"
    File.open(name + '_mteps_plt.txt', 'w') do |f|
        mteps = times.map{ |t| edges_count.to_f / t / 1e6 }
        mteps.each{ |p| f.puts p }
    end
end

impls = Set.new
graphs = Set.new
Dir.glob("*.raw") do |name|
    impl, graph = name.scan(/\w+/)
    impls.add(impl)
    graphs.add(graph)
    generate_data(impl, graph)
end
impls = impls.to_a

impls_names = impls.map{ |s| s = s.upcase; s = s.gsub('_', '-'); s}
colors = %w(red blue magenta orange green)
graphs.each do |graph|
    File.open(graph + '.plt', 'w') do |plot|
        plots = impls_names.size.times.map do |id| 
            filename = '"' + impls[id] + '-' + graph + '.txt' + '"'
            sp1 = 'with lines'
            color = 'lt rgb "' + colors[id] + '"'
            sp2 = 'lw 2'
            title = 'title "' + impls_names[id] + '"'
            [filename, sp1, color, sp2, title].join(' ')
        end.join(', ')

        plot.puts 'set terminal png'
        plot.puts "set output \"#{graph}.png\""
        plot.puts
        plot.puts 'set xtics ("1" 0, "2" 1, "4" 2, "8" 3, "16" 4)'
        plot.puts 'set xlabel "Threads count"'
        plot.puts 'set ylabel "Time, sec"'
        plot.puts "set title \"#{graph.upcase}\""
        plot.puts
        plot.puts 'plot ' + plots
    end
end

#p parse_times_from_file('errors')
#p best_times_from_directory("al-copy_rmat-22-32.raw")
