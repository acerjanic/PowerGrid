# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure(2) do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://atlas.hashicorp.com/search.
  config.vm.box = "boxcutter/ubuntu1504"

  config.vm.provider :aws do |aws, override|
    # config.vm.box = "dummy"
    aws.access_key_id = "#{ENV['AWS_ACCESS_KEY']}"
    aws.secret_access_key = "#{ENV['AWS_SECRET_KEY']}"
    aws.keypair_name = "VagrantKey"
    aws.ami = "emi-dc40a7ee"
    aws.instance_ready_timeout = 300
    aws.instance_type = "m1.medium"
    aws.tags = {
        "Name" => "VagrantPowerGrid",
    }
    aws.security_groups = ["Default"]
    aws.region = "eucalyptus"
    aws.endpoint = "http://bicloud.beckman.illinois.edu:8773/services/Eucalyptus"
    override.vm.box = "dummy"
    override.vm.box_url = "https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box"
    override.ssh.username = "ubuntu"
    override.ssh.private_key_path = "#{ENV['AWS_SSH_KEY']}"
    #override.vm.synced_folder ".", "/vagrant", type: "rsync"

    override.sync.host_folder = ""  #relative to the folder your Vagrantfile is in
    override.sync.guest_folder = "/vagrant" #relative to the vagrant home folder -> /home/vagrant
    #override.vm.synced_folder ".", "/vagrant", type: "nfs_guest"
    override.sshfs.enabled = false
    override.nfs.functional = false   
    override.sshfs.mount_on_guest = true
    override.sshfs.paths = { "" => "/vagrant" }
    override.sshfs.host_addr = "#{ENV['HOSTNAME']}"
    #override.ssh.private_key_path = "~/.ssh/id_rsa"
    #override.ssh.forward_agent = true
    override.ssh.insert_key = false
    override.ssh.private_key_path = ["#{ENV['AWS_SSH_KEY']}", '~/.ssh/id_rsa']
   end
  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  # config.vm.network "forwarded_port", guest: 80, host: 8080

  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
  # config.vm.network "private_network", ip: "192.168.33.10"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"

  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  # config.vm.provider "virtualbox" do |vb|
  #   # Display the VirtualBox GUI when booting the machine
  #   vb.gui = true
  #
  #   # Customize the amount of memory on the VM:
  #   vb.memory = "1024"
  # end

  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Define a Vagrant Push strategy for pushing to Atlas. Other push strategies
  # such as FTP and Heroku are also available. See the documentation at
  # https://docs.vagrantup.com/v2/push/atlas.html for more information.
  # config.push.define "atlas" do |push|
  #   push.app = "YOUR_ATLAS_USERNAME/YOUR_APPLICATION_NAME"
  # end

  # Enable provisioning with a shell script. Additional provisioners such as
  # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
  # documentation for more information about their specific syntax and use.
   config.vm.provision "shell", inline: <<-SHELL
     sudo apt-get update
  #  Installing ismrmrd
     sudo apt-get install -y libhdf5-serial-dev h5utils cmake cmake-curses-gui libboost-all-dev doxygen git libfftw3-dev g++
     git clone https://github.com/ismrmrd/ismrmrd
     cd ismrmrd/
     mkdir build
     cd build
     cmake ../
     make
     sudo make install
  #  Installing the PowerGrid Dependencies
     sudo apt-get install -y libmatio-dev libopenblas-dev libxerces-c-dev libarmadillo-dev xsdcxx

   SHELL
end

